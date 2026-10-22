# ============================================================================
# R/functions/water_access_functions.R
#
# Vulnerability 205: reliance on rainfall-dependent water sources (primary)
# and on unimproved sources (JMP definition, sensitivity).
#   prepare_unimproved_water_2022()  Census 2022 by police district
#     (DataFirst casd-2022-policedist v1) on the Stats SA police-district
#     polygons: rainfall-dependent (primary) and unimproved layers.
#   prepare_unimproved_water()       Census 2011 Community Profiles by SAL
#     (DataFirst 517, zaf-statssa-saccp-2011-v1) on SAL_SA_2013 in
#     mapdata-sa.zip (unimproved only; sensitivity).
# Needs sf, terra, haven.
# ============================================================================


# Census 2011 "main source of water" categories (source-of-water.dta)
#   wtrsrs_1  Regional/local water scheme      improved (piped)
#   wtrsrs_2  Borehole                         improved
#   wtrsrs_3  Spring                           UNIMPROVED (default)
#   wtrsrs_4  Rain water tank                  improved
#   wtrsrs_5  Dam/pool/stagnant water          UNIMPROVED
#   wtrsrs_6  River/stream                     UNIMPROVED
#   wtrsrs_7  Water vendor                     improved (delivered)
#   wtrsrs_8  Water tanker                     improved (delivered)
#   wtrsrs_9  Other                            excluded from the numerator
#   wtrsrs_10 Not applicable                   excluded (always 0)
#
# The census does not separate protected from unprotected springs. Most
# springs in the former Transkei are unprotected, so "spring" is counted as
# unimproved by default; set unimproved_cols = c("wtrsrs_5", "wtrsrs_6") for
# a strict surface-water-only version.
UNIMPROVED_WATER_COLS <- c("wtrsrs_3", "wtrsrs_5", "wtrsrs_6")

# Census 2011 "access to piped water" categories (piped-water.dta). Below the
# RDP basic-service standard = community stand > 200 m away, or no piped water.
BELOW_BASIC_PIPED_COLS <- c("pipwtr_4", "pipwtr_5", "pipwtr_6", "pipwtr_7")

# SANLC 2022 settlement classes: residential formal (47-50), residential
# informal (51-54), village (55-56), smallholdings (57-60).
SANLC_SETTLEMENT_CLASSES <- 47:60


# ----------------------------------------------------------------------------
# sal_household_counts(): read one Community Profile table
# ----------------------------------------------------------------------------
# Returns sal_code, sp_code, mp_code, lm_code (character) and the category
# columns, for SALs with a code (the file has one blank row).

sal_household_counts <- function(dta, prefix) {
  d <- haven::read_dta(dta)
  cols <- grep(paste0("^", prefix, "_[0-9]+$"), names(d), value = TRUE)
  if (!length(cols)) stop("No ", prefix, "_* columns in ", basename(dta))
  d <- as.data.frame(d[!is.na(d$sal_code),
                       c("sal_code", "sp_code", "mp_code", "mn_mdb_c", cols)])
  d$sal_code <- sprintf("%.0f", d$sal_code)
  d$sp_code  <- sprintf("%.0f", d$sp_code)
  d$mp_code  <- sprintf("%.0f", d$mp_code)
  names(d)[names(d) == "mn_mdb_c"] <- "lm_code"
  d
}


# ----------------------------------------------------------------------------
# pool_share(): share with small-area pooling
# ----------------------------------------------------------------------------
# A SAL with fewer than `min_hh` households takes the share of its sub-place
# (all SALs in that sub-place pooled); if the sub-place is still too small,
# its main place; then its local municipality. Pooling uses the national
# table, so sub-places cut by the ROI box are pooled in full.
#
# Stats SA randomly rounded Community Profile counts to base 3. Simulations
# of that rounding give a share error (SD) of ~0.15-0.2 at 15 households,
# ~0.08 at 30 and ~0.05 at 50; 30 is the default, 50 a stricter check.
#
# Returns a data frame: share, level ("sal", "sub_place", "main_place",
# "municipality"), pool_hh (households the share is based on). SALs with no
# households get NA.

pool_share <- function(num, den, groups, min_hh = 30) {
  share <- ifelse(den > 0, num / den, NA_real_)
  level <- ifelse(den > 0, "sal", NA_character_)
  pool  <- den

  for (g in names(groups)) {
    key <- groups[[g]]
    g_num <- ave(num, key, FUN = sum)
    g_den <- ave(den, key, FUN = sum)
    swap <- den > 0 & pool < min_hh & g_den > pool
    share[swap] <- g_num[swap] / g_den[swap]
    level[swap] <- g
    pool[swap]  <- g_den[swap]
  }
  data.frame(share = share, level = level, pool_hh = pool)
}


# ----------------------------------------------------------------------------
# prepare_unimproved_water(): SAL shares + household-weighted 1 km layers
# ----------------------------------------------------------------------------
#   mapdata_zip    : mapdata-sa.zip (only SAL_SA_2013 is extracted)
#   source_dta     : source-of-water.dta   (households by main water source)
#   piped_dta      : piped-water.dta       (households by piped-water access)
#   sanlc_file     : SANLC 2022 raster (the CRA extract)
#   roi            : the `roi` target (WGS84 bbox)
#   template       : the `template_grid` target (1 km, SA Albers)
#   out_gpkg       : SAL table with counts and shares (for maps and QA)
#   out_tif        : share of households using unimproved sources, 1 km
#   out_tif_piped  : share of households below the RDP piped-water standard,
#                    1 km (companion layer; optional indicator)
#   min_hh         : pooling threshold (households); see pool_share()
#   fine           : sub-cells per template cell side (10 -> 100 m)
#
# Method
#   1. Each SAL's share comes from pool_share() (SALs under `min_hh`
#      households borrow the share of their sub-place, main place or
#      municipality). The SAL keeps its own household count.
#   2. Households are spread onto 100 m sub-cells in proportion to the
#      fraction of each sub-cell under SANLC 2022 settlement classes
#      (47-60). SALs with households but no mapped settlement (scattered farm
#      dwellings) fall back to an even spread over the SAL; SALs too small to
#      contain a sub-cell centre go to the sub-cell at their centroid. No
#      households are lost.
#   3. Unimproved households (share x households) and total households are
#      summed to 1 km; the 1 km share is their ratio. Cells with no
#      households are NA.
#
# Settlement is mapped for 2022 and households counted in 2011, so new
# settlement since 2011 attracts some 2011 households; the SAL totals are
# unchanged.

prepare_unimproved_water <- function(mapdata_zip, source_dta, piped_dta,
                                     sanlc_file, roi, template,
                                     out_gpkg, out_tif, out_tif_piped = NULL,
                                     min_hh = 30, fine = 10,
                                     settlement_classes = SANLC_SETTLEMENT_CLASSES,
                                     unimproved_cols = UNIMPROVED_WATER_COLS,
                                     below_basic_cols = BELOW_BASIC_PIPED_COLS) {

  # 1. Household counts and pooled shares (national tables) -------------------
  src <- sal_household_counts(source_dta, "wtrsrs")
  src$hh_total      <- rowSums(src[grep("^wtrsrs_", names(src))])
  src$hh_unimproved <- rowSums(src[unimproved_cols])

  pip <- sal_household_counts(piped_dta, "pipwtr")
  pip$hh_piped_total <- rowSums(pip[paste0("pipwtr_", 1:7)])  # excl. unspecified/NA
  pip$hh_below_basic <- rowSums(pip[below_basic_cols])

  groups <- function(d) list(sub_place = d$sp_code, main_place = d$mp_code,
                             municipality = d$lm_code)

  p1 <- pool_share(src$hh_unimproved, src$hh_total, groups(src), min_hh)
  src$share_unimproved     <- p1$share
  src$pool_level_unimproved <- p1$level
  src$pool_hh_unimproved   <- p1$pool_hh

  p2 <- pool_share(pip$hh_below_basic, pip$hh_piped_total, groups(pip), min_hh)
  pip$share_below_basic     <- p2$share
  pip$pool_level_below_basic <- p2$level

  # 2. SAL polygons intersecting the ROI --------------------------------------
  exdir <- file.path(tempdir(), "census2011_sal")
  members <- paste0("mapdata-sa/SAL_SA_2013/SAL_SA_2013.",
                    c("shp", "shx", "dbf", "prj"))
  utils::unzip(mapdata_zip, files = members, exdir = exdir, junkpaths = TRUE)

  sal <- sf::st_read(
    file.path(exdir, "SAL_SA_2013.shp"), quiet = TRUE,
    wkt_filter = sf::st_as_text(sf::st_as_sfc(roi))
  ) |>
    sf::st_make_valid()
  sal$sal_code <- sprintf("%.0f", sal$SAL_CODE)
  sal <- sal[, c("sal_code", "SP_NAME", "MP_NAME", "MN_NAME", "DC_MDB_C",
                 "PR_NAME")]
  names(sal)[2:6] <- c("sp_name", "mp_name", "lm_name", "dm_code", "province")

  sal <- merge(sal, src[c("sal_code", "sp_code", "mp_code", "lm_code",
                          "hh_total", "hh_unimproved", "share_unimproved",
                          "pool_level_unimproved", "pool_hh_unimproved")],
               by = "sal_code", all.x = TRUE)
  sal <- merge(sal, pip[c("sal_code", "hh_piped_total", "hh_below_basic",
                          "share_below_basic", "pool_level_below_basic")],
               by = "sal_code", all.x = TRUE)
  if (anyNA(sal$hh_total)) {
    stop(sum(is.na(sal$hh_total)), " SALs in the ROI have no water-source table.")
  }

  # Households used for spreading: the SAL's own count x the (pooled) share
  sal$hh_unimproved_est <- ifelse(sal$hh_total > 0,
                                  sal$share_unimproved * sal$hh_total, 0)
  sal$hh_below_basic_est <- ifelse(sal$hh_piped_total > 0,
                                   sal$share_below_basic * sal$hh_piped_total, 0)

  dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(sal, out_gpkg, layer = "unimproved_water_sal2011",
               delete_dsn = TRUE, quiet = TRUE)

  # 3. Settlement weights on the 100 m sub-grid --------------------------------
  tmpl <- terra::rast(template)
  sub  <- terra::disagg(terra::rast(tmpl), fact = fine)   # empty 100 m grid

  sanlc  <- terra::rast(sanlc_file)
  settle <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
  settle <- terra::project(settle, sub, method = "average")  # fraction settled
  wt <- terra::values(settle, mat = FALSE)
  wt[is.na(wt)] <- 0

  v <- terra::project(terra::vect(sal), terra::crs(tmpl))
  v$zone <- seq_len(nrow(v))
  zone <- terra::values(terra::rasterize(v, sub, field = "zone"), mat = FALSE)

  in_zone    <- !is.na(zone)
  ncell_zone <- tabulate(zone[in_zone], nbins = nrow(v))
  wsum_zone  <- vapply(split(wt[in_zone], factor(zone[in_zone], seq_len(nrow(v)))),
                       sum, numeric(1))
  use_wt <- wsum_zone > 0                     # SAL has mapped settlement

  tiny <- which(ncell_zone == 0)
  tiny_cell <- if (length(tiny)) {
    terra::cellFromXY(sub, terra::crds(terra::centroids(v[tiny], inside = TRUE)))
  } else integer(0)

  # Share of each sub-cell in its SAL: settlement weight, or even spread
  z <- zone[in_zone]
  frac <- ifelse(use_wt[z], wt[in_zone] / wsum_zone[z], 1 / ncell_zone[z])

  spread <- function(counts) {
    x <- numeric(terra::ncell(sub))
    x[in_zone] <- counts[z] * frac
    for (i in seq_along(tiny)) {
      if (!is.na(tiny_cell[i])) x[tiny_cell[i]] <- x[tiny_cell[i]] + counts[tiny[i]]
    }
    terra::aggregate(terra::setValues(sub, x), fact = fine, fun = "sum")
  }

  to_share <- function(num, den, name) {
    r <- terra::ifel(den > 0, num / den, NA)
    r <- terra::mask(r, tmpl)
    names(r) <- name
    r
  }

  # 4. 1 km layers --------------------------------------------------------------
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  hh_tot  <- spread(v$hh_total)
  r_unimp <- to_share(spread(v$hh_unimproved_est), hh_tot,
                      "share_unimproved_water")
  dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r_unimp, out_tif, overwrite = TRUE, datatype = "FLT4S",
                     gdal = gdal)
  outs <- c(out_gpkg, out_tif)

  if (!is.null(out_tif_piped)) {
    r_pip <- to_share(spread(v$hh_below_basic_est), spread(v$hh_piped_total),
                      "share_below_basic_piped_water")
    terra::writeRaster(r_pip, out_tif_piped, overwrite = TRUE,
                       datatype = "FLT4S", gdal = gdal)
    outs <- c(outs, out_tif_piped)
  }

  # QA -------------------------------------------------------------------------
  pop <- v$hh_total > 0
  message(
    "SALs: ", nrow(v), " | with households: ", sum(pop),
    " | pooled (<", min_hh, " hh): ",
    paste(names(table(v$pool_level_unimproved[pop])),
          table(v$pool_level_unimproved[pop]), sep = " ", collapse = ", "),
    "\nPopulated SALs spread on settlement: ", sum(use_wt & pop),
    " | evenly (no mapped settlement): ", sum(!use_wt & pop & ncell_zone > 0),
    " | by centroid: ", sum(pop & ncell_zone == 0),
    "\nHouseholds: ", sum(v$hh_total),
    " | spread to grid: ", round(terra::global(hh_tot, "sum", na.rm = TRUE)[[1]]),
    " | ROI share unimproved: ",
    round(100 * sum(v$hh_unimproved_est) / sum(v$hh_total), 1), "%"
  )
  outs
}


# ============================================================================
# Census 2022: police districts (DataFirst)
# ============================================================================
# Source: Stats SA Census 2022 Community/Aggregate Statistics by police
# district (DataFirst, casd-2022-policedist v1), Stata tables
# c2022-watersource-policedist-v1.dta (qh06_*) and
# c2022-pipedwater-policedist-v1.dta (qh05_*); weighted household counts for
# all 1,167 districts. Police districts are the finest level these tables
# are published at.
#
# Boundaries: Stats SA police-district shapefile (PoliceDistrict.zip; SAPS
# component boundaries, 1,167 districts; fields PD_CODE, PD_Name, MN_MDB_C;
# see R/10_download/108_download_boundaries_police_districts.R). PD_Name
# matches the table's police_district one to one after police_key().
#
# Every household in a police district gets the district's share; households
# are placed on SANLC 2022 settlement so that 1 km cells straddling a
# district boundary are household-weighted.

# Census 2022 "Source of water" columns counted as unimproved: spring +
# dam/pool/stagnant water + river/stream (same definition as 2011)
UNIMPROVED_WATER_2022 <- c("qh06_spring", "qh06_dam_pool_stagnant", "qh06_flowing_river_stream")

# Census 2022 "Source of water" columns counted as rainfall-dependent: the
# unimproved sources + rain-water tank. These are the sources that fail when
# rainfall fails; water tankers (trucked in) are excluded.
RAINFALL_DEPENDENT_2022 <- c(UNIMPROVED_WATER_2022, "qh06_rainwater_tank")

# Census 2022 "Piped water" columns below the RDP standard (community stand
# more than 200 m away, or no piped water)
BELOW_BASIC_PIPED_2022 <- c("qh05_piped_500m", "qh05_piped_1000m",
                            "qh05_piped_morethan_1000m", "qh05_none")


# ----------------------------------------------------------------------------
# police_key(): name used for matching table rows to district polygons
# ----------------------------------------------------------------------------
# Upper case, "&" -> AND, other punctuation and repeated spaces removed:
#   "MIDDELBURG(EC)" -> "MIDDELBURG EC"; "Kwa Bhaca" -> "KWA BHACA"
police_key <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("&", " AND ", x, fixed = TRUE)
  x <- gsub("[^A-Z0-9]+", " ", x)
  trimws(gsub(" +", " ", x))
}


# ----------------------------------------------------------------------------
# read_district_table(): one DataFirst police-district table
# ----------------------------------------------------------------------------
# Returns district, district_key, local_muni and the numeric columns starting
# with `prefix` (e.g. "qh06_").

read_district_table <- function(dta, prefix) {
  d <- as.data.frame(haven::zap_labels(haven::read_dta(dta)))
  cols <- grep(paste0("^", prefix), names(d), value = TRUE)
  if (!length(cols)) stop(basename(dta), ": no columns starting with ", prefix)
  out <- data.frame(district = d$police_district, local_muni = d$local_muni,
                    stringsAsFactors = FALSE)
  for (cc in cols) out[[cc]] <- ifelse(is.na(d[[cc]]), 0, as.numeric(d[[cc]]))
  out$district_key <- police_key(out$district)
  if (anyDuplicated(out$district_key)) {
    stop(basename(dta), ": repeated districts: ",
         paste(unique(out$district[duplicated(out$district_key)]), collapse = ", "))
  }
  out
}


# ----------------------------------------------------------------------------
# read_police_districts(): Stats SA police-district polygons from the zip
# ----------------------------------------------------------------------------
# Reads the one .shp inside `police_zip` (no extraction), keeps districts
# intersecting `bbox_4326` (an sf bbox or sfc in EPSG:4326) and returns
# pd_code, pd_name, district_key, lm_code, then the geometry.

read_police_districts <- function(police_zip, bbox_4326 = NULL) {
  shp <- grep("\\.shp$", utils::unzip(police_zip, list = TRUE)$Name, value = TRUE)
  if (length(shp) != 1) stop(basename(police_zip), ": expected one .shp, found ",
                             length(shp))
  args <- list(dsn = file.path("/vsizip", normalizePath(police_zip, winslash = "/"), shp),
               quiet = TRUE)
  if (!is.null(bbox_4326)) {
    args$wkt_filter <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(bbox_4326)))
  }
  x <- do.call(sf::st_read, args) |> sf::st_make_valid()
  need <- c("PD_CODE", "PD_Name", "MN_MDB_C")
  if (!all(need %in% names(x))) {
    stop(basename(police_zip), ": missing fields ",
         paste(setdiff(need, names(x)), collapse = ", "))
  }
  out <- data.frame(pd_code = x$PD_CODE, pd_name = x$PD_Name,
                    district_key = police_key(x$PD_Name), lm_code = x$MN_MDB_C,
                    geometry = sf::st_geometry(x))
  out <- sf::st_as_sf(out)
  if (anyDuplicated(out$district_key)) {
    stop("Repeated police-district names: ",
         paste(unique(out$pd_name[duplicated(out$district_key)]), collapse = ", "))
  }
  out
}


# ----------------------------------------------------------------------------
# prepare_unimproved_water_2022(): police-district shares -> 1 km layers
# ----------------------------------------------------------------------------
#   source_dta    : DataFirst c2022-watersource-policedist-v1.dta
#   piped_dta     : DataFirst c2022-pipedwater-policedist-v1.dta, or NULL to
#                   skip the companion layer
#   police_zip    : Stats SA PoliceDistrict.zip
#   sanlc_file    : SANLC 2022 raster (the CRA extract)
#   template      : the `template_grid` target (1 km, SA Albers)
#   out_gpkg      : district table with counts and shares (maps and QA)
#   out_tif       : share of households using unimproved sources, 1 km
#                   (JMP definition; sensitivity layer)
#   out_tif_rain  : share of households using rainfall-dependent sources
#                   (unimproved + rain-water tank), 1 km (primary layer)
#   out_tif_piped : share below the RDP piped-water standard, 1 km
#   fine          : sub-cells per template cell side (10 -> 100 m)
#
# Method
#   1. Districts touching the template grid are joined to the tables by name;
#      every one must be in both tables (else stop).
#   2. District shares = unimproved (or rainfall-dependent) / all households
#      (weighted Census 2022 counts). Every household in the district gets
#      these shares.
#   3. District households are spread onto 100 m sub-cells by SANLC 2022
#      settlement (classes 47-60), evenly if a district has no mapped
#      settlement (spread_households(), crop_reliance_functions.R).
#   4. Numerator and total households are summed to 1 km; the 1 km share is
#      their ratio (cells with no households are NA). Inside one district the
#      value is constant; it changes only where cells straddle districts.

prepare_unimproved_water_2022 <- function(source_dta, piped_dta = NULL,
                                          police_zip, sanlc_file, template,
                                          out_gpkg, out_tif, out_tif_rain = NULL,
                                          out_tif_piped = NULL, fine = 10,
                                          settlement_classes = SANLC_SETTLEMENT_CLASSES,
                                          unimproved_cols = UNIMPROVED_WATER_2022,
                                          rain_cols = RAINFALL_DEPENDENT_2022,
                                          below_basic_cols = BELOW_BASIC_PIPED_2022) {

  tmpl <- terra::rast(template)

  # 1. Districts that touch the template grid --------------------------------
  grid_poly <- sf::st_as_sfc(sf::st_bbox(terra::ext(tmpl), crs = sf::st_crs(terra::crs(tmpl))))
  pd <- read_police_districts(police_zip, sf::st_transform(grid_poly, 4326))
  pd <- sf::st_transform(pd, terra::crs(tmpl))
  pd <- pd[lengths(sf::st_intersects(pd, grid_poly)) > 0, ]

  join <- function(tab, what) {
    miss <- setdiff(pd$district_key, tab$district_key)
    if (length(miss)) {
      stop("Police districts in the grid missing from ", what, ": ",
           paste(pd$pd_name[pd$district_key %in% miss], collapse = ", "))
    }
  }

  # 2. Source of water ---------------------------------------------------------
  src <- read_district_table(source_dta, "qh06_")
  stopifnot(all(c(unimproved_cols, rain_cols) %in% names(src)))
  src$hh_total      <- rowSums(src[grep("^qh06_", names(src))])
  src$hh_unimproved <- rowSums(src[unimproved_cols])
  src$hh_rain_dep   <- rowSums(src[rain_cols])
  join(src, basename(source_dta))
  pd <- merge(pd, src[c("district_key", "local_muni", "hh_total", "hh_unimproved", "hh_rain_dep")],
              by = "district_key", all.x = TRUE)
  pd$share_unimproved <- ifelse(pd$hh_total > 0, pd$hh_unimproved / pd$hh_total, NA)
  pd$share_rain_dep   <- ifelse(pd$hh_total > 0, pd$hh_rain_dep / pd$hh_total, NA)

  # 3. Piped water (companion) -------------------------------------------------
  if (!is.null(piped_dta)) {
    pip <- read_district_table(piped_dta, "qh05_")
    stopifnot(all(below_basic_cols %in% names(pip)))
    pip$hh_piped_total <- rowSums(pip[grep("^qh05_", names(pip))])
    pip$hh_below_basic <- rowSums(pip[below_basic_cols])
    join(pip, basename(piped_dta))
    pd <- merge(pd, pip[c("district_key", "hh_piped_total", "hh_below_basic")],
                by = "district_key", all.x = TRUE)
    pd$share_below_basic <- ifelse(pd$hh_piped_total > 0,
                                   pd$hh_below_basic / pd$hh_piped_total, NA)
  }

  pd <- pd[, c("pd_code", "pd_name", "district_key", "lm_code", "local_muni",
               setdiff(names(pd), c("pd_code", "pd_name", "district_key",
                                    "lm_code", "local_muni")))]
  dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(pd, out_gpkg, layer = "unimproved_water_police2022",
               delete_dsn = TRUE, quiet = TRUE)

  # 4. Spread households on settlement and sum to 1 km -------------------------
  sub <- terra::disagg(terra::rast(tmpl), fact = fine)
  sanlc  <- terra::rast(sanlc_file)
  settle <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
  settle <- terra::project(settle, sub, method = "average")

  v <- terra::vect(pd)
  for (f in intersect(c("hh_total", "hh_unimproved", "hh_rain_dep", "hh_piped_total",
                        "hh_below_basic"), names(v))) {
    x <- v[[f, drop = TRUE]]; x[is.na(x)] <- 0; v[[f]] <- x
  }
  spread_1km <- function(field) {
    sp <- spread_households(v, sub, settle, field)
    list(r = terra::aggregate(sp$surface, fact = fine, fun = "sum"), n = sp$n)
  }
  to_share <- function(num, den, name) {
    r <- terra::ifel(den > 0, num / den, NA)
    r <- terra::mask(r, tmpl)
    names(r) <- name
    r
  }

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  tot <- spread_1km("hh_total")
  r_unimp <- to_share(spread_1km("hh_unimproved")$r, tot$r, "share_unimproved_water_2022")
  dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r_unimp, out_tif, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  outs <- c(out_gpkg, out_tif)

  if (!is.null(out_tif_rain)) {
    r_rain <- to_share(spread_1km("hh_rain_dep")$r, tot$r, "share_rainfall_dependent_water_2022")
    dir.create(dirname(out_tif_rain), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(r_rain, out_tif_rain, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
    outs <- c(outs, out_tif_rain)
  }

  if (!is.null(piped_dta) && !is.null(out_tif_piped)) {
    r_pip <- to_share(spread_1km("hh_below_basic")$r, spread_1km("hh_piped_total")$r,
                      "share_below_basic_piped_water_2022")
    terra::writeRaster(r_pip, out_tif_piped, overwrite = TRUE, datatype = "FLT4S",
                       gdal = gdal)
    outs <- c(outs, out_tif_piped)
  }

  # QA ---------------------------------------------------------------------------
  in_grid <- terra::global(terra::mask(tot$r, tmpl), "sum", na.rm = TRUE)[[1]]
  message(
    "Police districts touching the grid: ", nrow(pd),
    " | spread on settlement: ", tot$n[["settlement"]],
    ", evenly: ", tot$n[["even"]], ", by centroid: ", tot$n[["centroid"]],
    "\nHouseholds in these districts: ", round(sum(pd$hh_total)),
    " | inside the CRA extent: ", round(in_grid),
    "\nAll households in these districts: unimproved ",
    round(100 * sum(pd$hh_unimproved) / sum(pd$hh_total), 1), "%, rainfall-dependent ",
    round(100 * sum(pd$hh_rain_dep) / sum(pd$hh_total), 1), "%"
  )
  outs
}
