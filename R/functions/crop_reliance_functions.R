# ============================================================================
# R/functions/crop_reliance_functions.R
#
# Vulnerability indicator: reliance on rainfed crops, measured as SANLC 2022
# subsistence / small-scale cropland (class 41) per Census 2011 household.
#
# Two versions are produced from the same household and cropland surfaces:
#
#   1. Neighbourhood ratio on the 1 km template grid (the indicator):
#      cropland within `radius_m` of the cell / households within `radius_m`.
#      Households live in villages and farm fields around them, often in a
#      different SAL, so a per-SAL ratio mostly reflects where boundaries
#      fall. A radius of 3-5 km (walking distance to fields) matches
#      households to the fields they farm.
#   2. 2011 wards (companion layer and check): cropland and households
#      summed per ward. Stats SA's own 2011 wards (WD_SA_2011, same zip as
#      the SALs) are used; they are identical to the MDB 2011 wards.
#
# Why not per SAL: in a test (29 Sep 2026) 61% of households lived in SALs
# with < 0.01 ha/household, while a few near-empty SALs had > 80 ha. Also,
# SAL_SA_2013 covers only ~66% of the study area; the gaps are almost
# unpopulated (2% of settlement) but hold ~59% of class 41 cropland. Cropland
# is therefore always taken from the raster over the whole grid or ward,
# never summed over SAL polygons.
#
# Needs: sf, terra, haven. Uses sal_household_counts() and
# SANLC_SETTLEMENT_CLASSES from water_access_functions.R (tar_source() loads
# both files; the stage script sources both).
# ============================================================================


# SANLC 2022 class 41 = subsistence / small-scale annual crops. In this area
# these are rainfed. Add 40 (commercial annual crops, rain-fed/dryland) for a
# version that includes commercial dryland farming.
SANLC_SUBSISTENCE_CROP_CLASSES <- 41


# ----------------------------------------------------------------------------
# read_census2011_layer(): one layer from mapdata-sa.zip, clipped to the ROI
# ----------------------------------------------------------------------------
#   layer : e.g. "SAL_SA_2013" (small areas) or "WD_SA_2011" (wards)

read_census2011_layer <- function(mapdata_zip, layer, roi) {
  exdir <- file.path(tempdir(), paste0("census2011_", layer))
  members <- paste0("mapdata-sa/", layer, "/", layer, ".",
                    c("shp", "shx", "dbf", "prj"))
  utils::unzip(mapdata_zip, files = members, exdir = exdir, junkpaths = TRUE)
  sf::st_read(
    file.path(exdir, paste0(layer, ".shp")), quiet = TRUE,
    wkt_filter = sf::st_as_text(sf::st_as_sfc(roi))
  ) |>
    sf::st_make_valid()
}


# ----------------------------------------------------------------------------
# spread_households(): SAL household counts -> sub-grid surface
# ----------------------------------------------------------------------------
# Same method as prepare_unimproved_water(): households are spread onto the
# sub-grid in proportion to SANLC settlement cover; SALs with households but
# no mapped settlement are spread evenly; SALs too small to hold a sub-cell
# centre go to the sub-cell at their centroid. No households are lost.
#
#   v      : SAL SpatVector (template CRS) with a household column
#   sub    : empty sub-grid (e.g. 100 m)
#   settle : settlement fraction on `sub` (0-1, NA treated as 0)
#   field  : household column in `v`
#
# Returns a list: surface (SpatRaster on `sub`) and n (counts of SALs per
# spreading method, for QA).

spread_households <- function(v, sub, settle, field = "hh_total") {
  wt <- terra::values(settle, mat = FALSE)
  wt[is.na(wt)] <- 0

  v$zone <- seq_len(nrow(v))
  zone <- terra::values(terra::rasterize(v, sub, field = "zone"), mat = FALSE)

  in_zone    <- !is.na(zone)
  ncell_zone <- tabulate(zone[in_zone], nbins = nrow(v))
  wsum_zone  <- vapply(split(wt[in_zone], factor(zone[in_zone], seq_len(nrow(v)))),
                       sum, numeric(1))
  use_wt <- wsum_zone > 0

  tiny <- which(ncell_zone == 0)
  tiny_cell <- if (length(tiny)) {
    terra::cellFromXY(sub, terra::crds(terra::centroids(v[tiny], inside = TRUE)))
  } else integer(0)

  counts <- v[[field, drop = TRUE]]
  counts[is.na(counts)] <- 0

  z <- zone[in_zone]
  frac <- ifelse(use_wt[z], wt[in_zone] / wsum_zone[z], 1 / ncell_zone[z])

  x <- numeric(terra::ncell(sub))
  x[in_zone] <- counts[z] * frac
  for (i in seq_along(tiny)) {
    if (!is.na(tiny_cell[i])) x[tiny_cell[i]] <- x[tiny_cell[i]] + counts[tiny[i]]
  }

  pop <- counts > 0
  list(
    surface = terra::setValues(sub, x),
    n = c(settlement = sum(use_wt & pop),
          even       = sum(!use_wt & pop & ncell_zone > 0),
          centroid   = sum(pop & ncell_zone == 0),
          lost       = sum(pop & ncell_zone == 0 & is.na(tiny_cell[match(seq_along(pop), tiny)])))
  )
}


# ----------------------------------------------------------------------------
# prepare_cropland_per_household()
# ----------------------------------------------------------------------------
#   mapdata_zip    : mapdata-sa.zip (SAL_SA_2013 and WD_SA_2011 are read)
#   households_dta : any Census 2011 Community Profile household table by
#                    SAL; households = sum of its categories. Default use:
#                    source-of-water.dta ("wtrsrs" columns).
#   hh_prefix      : column prefix of that table
#   sanlc_file     : SANLC 2022 raster (the CRA extract)
#   roi            : the `roi` target (WGS84 bbox)
#   template       : the `template_grid` target (1 km, SA Albers)
#   out_tif        : neighbourhood ratio, ha/household, 1 km (the indicator)
#   out_gpkg_wards : ward table: households, cropland, ratio
#   out_tif_wards  : ward ratio on the 1 km grid (optional companion layer)
#   crop_classes   : SANLC classes counted as rainfed smallholder cropland
#   radius_m       : neighbourhood radius in metres (sensitivity: 3000-5000)
#   min_hh_nbhd    : neighbourhoods with fewer households than this are NA.
#                    50 keeps the extreme tail down (max ~32 ha/household
#                    at 50 vs ~46 at 10; household-weighted p10/p50/p90
#                    unchanged at 0.07/0.75/1.65 for r = 5 km)
#   mask_to_households : TRUE = only cells with households get a value
#   fine           : sub-cells per template cell side (10 -> 100 m)
#
# Edge effect: cropland and households outside the template grid are not
# counted, so cells within `radius_m` of the grid edge see only part of
# their neighbourhood. Both terms are cut in the same way, so the ratio is
# less affected than either total; the grid edge lies outside the study
# area in most places.
#
# Time mismatch: 2011 households, 2022 cropland and settlement.

prepare_cropland_per_household <- function(mapdata_zip, households_dta,
                                           sanlc_file, roi, template,
                                           out_tif, out_gpkg_wards,
                                           out_tif_wards = NULL,
                                           hh_prefix = "wtrsrs",
                                           crop_classes = SANLC_SUBSISTENCE_CROP_CLASSES,
                                           radius_m = 5000,
                                           min_hh_nbhd = 50,
                                           mask_to_households = TRUE,
                                           fine = 10,
                                           settlement_classes = SANLC_SETTLEMENT_CLASSES) {

  tmpl <- terra::rast(template)
  sub  <- terra::disagg(terra::rast(tmpl), fact = fine)     # empty sub-grid
  sub_ha <- prod(terra::res(sub)) / 1e4                     # 1 ha at 100 m

  # 1. Households per SAL (full-count Community Profile table) --------------
  hh <- sal_household_counts(households_dta, hh_prefix)
  hh$hh_total <- rowSums(hh[grep(paste0("^", hh_prefix, "_"), names(hh))])

  sal <- read_census2011_layer(mapdata_zip, "SAL_SA_2013", roi)
  sal$sal_code <- sprintf("%.0f", sal$SAL_CODE)
  sal <- merge(sal[, "sal_code"], hh[c("sal_code", "hh_total")],
               by = "sal_code", all.x = TRUE)
  if (anyNA(sal$hh_total)) {
    stop(sum(is.na(sal$hh_total)), " SALs in the ROI have no household count.")
  }

  # 2. Settlement and cropland on the sub-grid -------------------------------
  sanlc  <- terra::rast(sanlc_file)
  settle <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
  settle <- terra::project(settle, sub, method = "average")

  crop <- terra::ifel(is.na(sanlc), NA, sanlc %in% crop_classes) * 1
  crop <- terra::project(crop, sub, method = "average")     # fraction cropped
  crop <- terra::ifel(is.na(crop), 0, crop) * sub_ha        # ha per sub-cell

  v <- terra::project(terra::vect(sal), terra::crs(tmpl))
  hh_spread <- spread_households(v, sub, settle, "hh_total")
  hh_sub <- hh_spread$surface

  # 3. Neighbourhood ratio on the 1 km grid -----------------------------------
  crop_1km <- terra::aggregate(crop,   fact = fine, fun = "sum")
  hh_1km   <- terra::aggregate(hh_sub, fact = fine, fun = "sum")

  w <- terra::focalMat(tmpl, radius_m, type = "circle")
  w[w > 0] <- 1
  crop_nb <- terra::focal(crop_1km, w, fun = "sum", na.rm = TRUE)
  hh_nb   <- terra::focal(hh_1km,   w, fun = "sum", na.rm = TRUE)

  ratio <- terra::ifel(hh_nb >= min_hh_nbhd, crop_nb / hh_nb, NA)
  if (mask_to_households) ratio <- terra::ifel(hh_1km > 0, ratio, NA)
  ratio <- terra::mask(ratio, tmpl)
  names(ratio) <- sprintf("cropland_ha_per_household_r%gkm", radius_m / 1000)

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(ratio, out_tif, overwrite = TRUE, datatype = "FLT4S",
                     gdal = gdal)
  outs <- out_tif

  # 4. 2011 wards --------------------------------------------------------------
  # Summing the sub-grid surfaces by ward allocates SAL households that
  # straddle a ward boundary by settlement, not by largest overlap.
  wards <- read_census2011_layer(mapdata_zip, "WD_SA_2011", roi)
  wards <- sf::st_transform(wards, terra::crs(tmpl))
  wards <- data.frame(
    ward_id  = as.character(wards$WARD_ID),
    ward_no  = wards$WARDNO,
    lm_code  = wards$CAT_B,
    province = wards$PROVINCE,
    geometry = sf::st_geometry(wards)
  ) |> sf::st_as_sf()

  wv <- terra::vect(wards)
  wv$zone <- seq_len(nrow(wv))
  wzone <- terra::rasterize(wv, sub, field = "zone")
  z_hh   <- terra::zonal(hh_sub, wzone, fun = "sum", na.rm = TRUE)
  z_crop <- terra::zonal(crop,   wzone, fun = "sum", na.rm = TRUE)

  wards$households <- 0
  wards$households[z_hh[[1]]] <- z_hh[[2]]
  wards$cropland_ha <- 0
  wards$cropland_ha[z_crop[[1]]] <- z_crop[[2]]
  wards$ha_per_household <- ifelse(wards$households > 0,
                                   wards$cropland_ha / wards$households, NA)
  wards$area_km2 <- as.numeric(sf::st_area(wards)) / 1e6

  # Share of each ward inside the grid and inside the study area. Wards cut
  # by the grid edge have partial totals.
  grid_poly  <- sf::st_as_sfc(sf::st_bbox(terra::ext(tmpl), crs = sf::st_crs(wards)))
  study_poly <- sf::st_as_sf(terra::as.polygons(tmpl)) |> sf::st_union()
  frac_in <- function(p) {
    a <- sf::st_area(sf::st_intersection(sf::st_geometry(wards), p))
    ids <- sf::st_intersects(wards, p, sparse = FALSE)[, 1]
    out <- rep(0, nrow(wards))
    out[ids] <- as.numeric(a) / (wards$area_km2[ids] * 1e6)
    pmin(out, 1)
  }
  wards$share_in_grid  <- frac_in(grid_poly)
  wards$share_in_study <- frac_in(study_poly)

  dir.create(dirname(out_gpkg_wards), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(wards, out_gpkg_wards, layer = "cropland_per_household_wards2011",
               delete_dsn = TRUE, quiet = TRUE)
  outs <- c(outs, out_gpkg_wards)

  if (!is.null(out_tif_wards)) {
    wr <- terra::rasterize(terra::vect(wards), tmpl, field = "ha_per_household")
    wr <- terra::mask(wr, tmpl)
    names(wr) <- "cropland_ha_per_household_ward2011"
    terra::writeRaster(wr, out_tif_wards, overwrite = TRUE, datatype = "FLT4S",
                       gdal = gdal)
    outs <- c(outs, out_tif_wards)
  }

  # QA -------------------------------------------------------------------------
  in_study <- !is.na(terra::values(tmpl, mat = FALSE))
  rv <- terra::values(ratio, mat = FALSE)[in_study]
  hv <- terra::values(hh_1km, mat = FALSE)[in_study]
  ok <- is.finite(rv)
  o  <- order(rv[ok])
  cw <- cumsum(hv[ok][o]) / sum(hv[ok])
  wq <- function(p) rv[ok][o][which(cw >= p)[1]]

  sw <- wards$share_in_study > 0.5
  message(
    "SALs: ", nrow(v), " | households spread by settlement: ",
    hh_spread$n[["settlement"]], ", evenly: ", hh_spread$n[["even"]],
    ", by centroid: ", hh_spread$n[["centroid"]],
    "\nHouseholds (SALs in ROI): ", sum(v$hh_total),
    " | on grid: ", round(sum(terra::values(hh_1km), na.rm = TRUE)),
    " | cropland on grid: ", round(sum(terra::values(crop_1km), na.rm = TRUE)), " ha",
    "\nNeighbourhood ratio (r = ", radius_m / 1000, " km), household-weighted",
    " p10/p50/p90: ", paste(round(c(wq(0.1), wq(0.5), wq(0.9)), 2), collapse = " / "),
    " ha/household; study-area cells with a value: ",
    round(100 * mean(ok), 1), "%",
    "\nWards > 50% in study area: ", sum(sw), " | households: ",
    round(sum(wards$households[sw])), " | cropland: ",
    round(sum(wards$cropland_ha[sw])), " ha | overall: ",
    round(sum(wards$cropland_ha[sw]) / sum(wards$households[sw]), 2),
    " ha/household"
  )
  outs
}
