# ============================================================================
# R/functions/settlement_functions.R
#
# Exposure indicator 211 "Settlements": rural households per km2, from the
# DWS settlement layer (Water Services Knowledge System, WSKS).
#
# The DWS layer has one polygon per settlement (median ~0.8 km2) with a
# household count for the financial year (FinYear "April 2026"), modelled by
# DWS from Census 2011 and updated annually. Checked 9 Oct 2026 against
# Census 2022 police-district totals: r = 0.95, overall ratio 0.99. DWS
# misses ~10% of Census 2022 households, mostly in towns and on commercial
# farms (e.g. Ngangelizwe, Mthatha, Ugie, Richmond-KZN), i.e. mostly outside
# the target smallholder communities.
#
# Why households and not SANLC settlement area: the two agree on WHERE
# people live (rank correlation 0.91 per 1 km cell), but households per ha
# of SANLC settlement range 0.6-8.9 between cells (p10-p90), and SANLC maps
# ~90% of settlement here as "residential formal", so it cannot separate
# villages from towns. DWS has a Rural/Urban flag (Category).
#
# Do not use the DWS population, growth or "PopDensity" fields (PopDensity is
# persons per household, and some LMs have 6.5-8 persons per household).
# The water-source fields are one profile per local municipality, so they
# are not used either (see claude/dws-settlements-water-access.md).
#
# Needs: sf, terra. Uses spread_households() (crop_reliance_functions.R) and
# SANLC_SETTLEMENT_CLASSES (water_access_functions.R); tar_source() loads
# all three files, the stage script sources them.
# ============================================================================


# Fields the functions below rely on
DWS_SETTLEMENT_FIELDS <- c("UserNr", "CommunityName", "Category", "LMCode",
                           "LMName", "TotalHouseholds")


# ----------------------------------------------------------------------------
# read_dws_settlements(): the DWS settlement polygons
# ----------------------------------------------------------------------------
#   dws_gdb : the DWS file geodatabase (data/raw/dws/Umzimvubu_Settlements.gdb)
#   layer   : layer name; NULL takes the only layer in the gdb
#   crs     : CRS to transform to (NULL keeps EPSG:4326)
#
# Geometries are repaired with st_make_valid() (4 of 5,739 are invalid) and
# empty ones dropped. Some polygons overlap slightly (~180 pairs, mostly
# slivers); spread_households() assigns an overlapping sub-cell to one of
# them, and no households are lost.

read_dws_settlements <- function(dws_gdb, layer = NULL, crs = NULL) {
  layers <- sf::st_layers(dws_gdb)$name
  if (is.null(layer)) {
    if (length(layers) != 1) {
      stop(basename(dws_gdb), " has ", length(layers), " layers (",
           paste(layers, collapse = ", "), "); set `layer`.")
    }
    layer <- layers
  }
  x <- sf::st_read(dws_gdb, layer = layer, quiet = TRUE)

  miss <- setdiff(DWS_SETTLEMENT_FIELDS, names(x))
  if (length(miss)) {
    stop(basename(dws_gdb), ": missing fields ", paste(miss, collapse = ", "))
  }

  x <- sf::st_make_valid(x)
  x <- x[!sf::st_is_empty(x), ]
  x$TotalHouseholds[is.na(x$TotalHouseholds)] <- 0
  if (!is.null(crs)) x <- sf::st_transform(x, crs)
  x
}


# ----------------------------------------------------------------------------
# prepare_settlement_households(): DWS households -> households per km2, 1 km
# ----------------------------------------------------------------------------
#   dws_gdb     : DWS settlement geodatabase
#   sanlc_file  : SANLC 2022 raster (the CRA extract)
#   template    : the `template_grid` target (1 km, SA Albers)
#   out_tif     : the indicator: households per km2 in settlements of
#                 `categories` (default rural only), 1 km
#   out_tif_all : optional companion: all households (rural + urban) per km2
#   out_gpkg    : optional settlement table (polygons touching the study
#                 area) for maps and QA
#   categories  : DWS Category values in the indicator ("Rural", "Urban")
#   fine        : sub-cells per template cell side (10 -> 100 m)
#   use_sanlc   : TRUE places each settlement's households on the SANLC
#                 settlement pixels inside its polygon; FALSE spreads them
#                 evenly over the polygon
#   pad_cells   : template cells added around the grid while spreading, so
#                 settlements crossing the grid edge keep the households
#                 that belong outside it (10 -> 10 km)
#
# Method
#   1. Settlement polygons touching the (padded) template grid are read and
#      put in the template CRS.
#   2. Each settlement's households are spread onto a 100 m sub-grid with
#      spread_households(): in proportion to SANLC 2022 settlement cover
#      (classes 47-60) inside the polygon, evenly if the polygon has no
#      mapped settlement, and to the sub-cell at its centroid if it is too
#      small to hold a sub-cell centre. Totals are kept exactly.
#   3. Sub-cells are summed to 1 km and divided by the cell area (km2).
#      Template cells without households are 0 (no households in a DWS
#      settlement), not NA; cells outside the CRA extent are NA.
#
# Returns the files written, indicator first.

prepare_settlement_households <- function(dws_gdb, sanlc_file, template,
                                          out_tif, out_tif_all = NULL,
                                          out_gpkg = NULL,
                                          categories = "Rural",
                                          fine = 10, use_sanlc = TRUE,
                                          pad_cells = 10,
                                          settlement_classes = SANLC_SETTLEMENT_CLASSES,
                                          layer = NULL) {

  tmpl <- terra::rast(template)
  tcrs <- terra::crs(tmpl)
  padded <- terra::extend(terra::rast(tmpl), pad_cells)

  # 1. Settlements touching the padded grid ----------------------------------
  grid_poly <- sf::st_as_sfc(sf::st_bbox(terra::ext(padded), crs = sf::st_crs(tcrs)))
  s <- read_dws_settlements(dws_gdb, layer = layer, crs = tcrs)
  s <- s[lengths(sf::st_intersects(s, grid_poly)) > 0, ]

  bad <- setdiff(categories, unique(s$Category))
  if (length(bad)) {
    stop("Category value(s) not in the DWS layer: ", paste(bad, collapse = ", "),
         ". Values present: ", paste(sort(unique(s$Category)), collapse = ", "))
  }

  # 2. Sub-grid and settlement weights ---------------------------------------
  sub <- terra::disagg(padded, fact = fine)
  if (use_sanlc) {
    sanlc  <- terra::rast(sanlc_file)
    settle <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
    settle <- terra::project(settle, sub, method = "average")
  } else {
    settle <- terra::init(sub, 0)    # no weights: every polygon spread evenly
  }

  cell_km2 <- prod(terra::res(tmpl)) / 1e6

  per_km2 <- function(sel, name) {
    v  <- terra::vect(sel)
    sp <- spread_households(v, sub, settle, field = "TotalHouseholds")
    r  <- terra::aggregate(sp$surface, fact = fine, fun = "sum") / cell_km2
    r  <- terra::crop(r, tmpl)
    r  <- terra::ifel(is.na(tmpl), NA, r)
    names(r) <- name
    list(r = r, n = sp$n)
  }

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  rural <- s[s$Category %in% categories, ]
  ind <- per_km2(rural, paste0(tolower(paste(categories, collapse = "_")),
                               "_households_per_km2"))
  dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(ind$r, out_tif, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  outs <- out_tif

  if (!is.null(out_tif_all)) {
    all_hh <- per_km2(s, "all_households_per_km2")
    dir.create(dirname(out_tif_all), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(all_hh$r, out_tif_all, overwrite = TRUE, datatype = "FLT4S",
                       gdal = gdal)
    outs <- c(outs, out_tif_all)
  }

  # 3. Settlement table (polygons touching the CRA extent) -------------------
  study <- terra::as.polygons(tmpl) |> sf::st_as_sf() |> sf::st_union()
  in_study <- lengths(sf::st_intersects(s, study)) > 0
  if (!is.null(out_gpkg)) {
    keep <- intersect(c("FinYear", "Province", "DMCode", "DMName", "LMCode", "LMName",
                        "UserNr", "CommunityName", "Category", "TotalHouseholds"),
                      names(s))
    g <- s[in_study, keep]
    g$area_km2 <- as.numeric(sf::st_area(g)) / 1e6
    g$households_per_km2 <- ifelse(g$area_km2 > 0, g$TotalHouseholds / g$area_km2, NA)
    g$in_indicator <- g$Category %in% categories
    dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
    sf::st_write(g, out_gpkg, layer = "dws_settlements", delete_dsn = TRUE, quiet = TRUE)
    outs <- c(outs, out_gpkg)
  }

  # QA ---------------------------------------------------------------------------
  vals <- terra::values(terra::mask(ind$r, tmpl), mat = FALSE)
  vals <- vals[is.finite(vals)]
  pos  <- vals[vals > 0]
  q <- if (length(pos)) round(stats::quantile(pos, c(0.5, 0.9, 0.99, 1)), 1) else NA
  # Check: households of indicator settlements with a point inside the study
  # area. Should be close to the gridded total (polygons crossing the
  # boundary make the small difference).
  pts <- suppressWarnings(sf::st_point_on_surface(sf::st_geometry(rural)))
  expected <- sum(rural$TotalHouseholds[lengths(sf::st_intersects(pts, study)) > 0])
  message(
    "DWS settlements near the grid: ", nrow(s), " (", sum(in_study),
    " touch the CRA extent); in the indicator (", paste(categories, collapse = ", "),
    "): ", nrow(rural),
    "\nSpread on SANLC settlement: ", ind$n[["settlement"]], ", evenly: ", ind$n[["even"]],
    ", by centroid: ", ind$n[["centroid"]],
    ", centre outside the padded grid (dropped): ", ind$n[["lost"]],
    "\nHouseholds in the indicator inside the CRA extent: ",
    round(sum(vals) * cell_km2), " (settlements with their centre inside: ",
    round(expected), ")",
    "\nTemplate cells with households: ", length(pos), " of ", length(vals),
    " (", round(100 * length(pos) / length(vals)), "%)",
    "\nHouseholds per km2 in those cells: median ", q[1], ", p90 ", q[2],
    ", p99 ", q[3], ", max ", q[4]
  )
  outs
}
