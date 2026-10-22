# ============================================================================
# R/functions/income_functions.R
#
# Vulnerability 210 "Lower-income areas": share of households with an annual
# income <= R19,600 (Census 2011 by SAL, DataFirst 517,
# annual-household-income.dta), spread on SANLC 2022 settlement and summed
# to the 1 km template. Template cells without households are NA (filled in
# stage 30 through the register's na_fill).
#
# prepare_low_income_2011() writes the 1 km share and a SAL table. Needs sf,
# terra, haven. Uses sal_household_counts(), pool_share(),
# SANLC_SETTLEMENT_CLASSES (water_access_functions.R) and spread_households()
# (crop_reliance_functions.R).
# ============================================================================


# Census 2011 annual household income (inchh_*):
#   1 No income; 2 R1-4,800; 3 R4,801-9,600; 4 R9,601-19,600;
#   5 R19,601-38,200; 6 R38,201-76,400; 7 R76,401-153,800; ...;
#   12 R2,457,601+; 13 Unspecified
# <= R19,600/yr (main) is close to the 2011 food poverty line for a
# four-person household; <= R38,200 (sensitivity) is close to the upper-bound
# poverty line but leaves little contrast here (84-88% in the rural LMs).
LOW_INCOME_COLS    <- paste0("inchh_", 1:4)
LOW_INCOME_COLS_UB <- paste0("inchh_", 1:5)
NO_INCOME_COLS     <- "inchh_1"
INCOME_UNSPEC_COL  <- "inchh_13"

# ----------------------------------------------------------------------------
# grid_setup(): template, padded grid, 100 m sub-grid, grid bbox in WGS84
# ----------------------------------------------------------------------------
grid_setup <- function(template, pad_cells, fine) {
  tmpl <- terra::rast(template)
  pad  <- terra::extend(terra::rast(tmpl), pad_cells)
  bb   <- sf::st_as_sfc(sf::st_bbox(terra::ext(pad), crs = sf::st_crs(terra::crs(tmpl))))
  list(tmpl = tmpl, pad = pad, sub = terra::disagg(pad, fact = fine),
       bb = bb, bb_ll = sf::st_transform(bb, 4326),
       km2 = prod(terra::res(tmpl)) / 1e6)
}

settlement_weights <- function(sanlc_file, sub, classes) {
  s <- terra::rast(sanlc_file)
  s <- terra::ifel(is.na(s), NA, s %in% classes) * 1
  terra::project(s, sub, method = "average")
}

to_template <- function(r, tmpl, name) {
  r <- terra::crop(r, tmpl)
  r <- terra::mask(r, tmpl)
  names(r) <- name
  r
}

write_1km <- function(r, f) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r, f, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  f
}


# ----------------------------------------------------------------------------
# prepare_low_income_2011(): Census 2011 SAL income -> 1 km layers
# ----------------------------------------------------------------------------
#   mapdata_zip  : mapdata-sa.zip (SAL_SA_2013)
#   income_dta   : annual-household-income.dta
#   sanlc_file   : SANLC 2022 raster (CRA extract)
#   template     : template grid (1 km, SA Albers)
#   out_share    : low-income share of households (indicator)
#   out_gpkg     : SAL table (counts, shares at three cut-offs, pooling level)
#   low_cols     : income columns counted as low income
#   min_hh       : pooling threshold (see pool_share())
#
# Method
#   1. Share = low-income / (all - unspecified) per SAL, pooled to sub-place,
#      main place, then LM for SALs < min_hh households (national table).
#   2. Low-income households = pooled share x the SAL's own household count.
#   3. SAL households spread on SANLC 2022 settlement (classes 47-60) at
#      100 m on a grid padded by pad_cells (edge SALs keep their outside
#      households), evenly if a SAL has no mapped settlement; summed to 1 km.

prepare_low_income_2011 <- function(mapdata_zip, income_dta, sanlc_file, template,
                                    out_share, out_gpkg = NULL,
                                    low_cols = LOW_INCOME_COLS, min_hh = 30,
                                    fine = 10, pad_cells = 10,
                                    settlement_classes = SANLC_SETTLEMENT_CLASSES) {

  g <- grid_setup(template, pad_cells, fine)

  # 1. Counts and pooled shares (national table) -------------------------------
  d  <- sal_household_counts(income_dta, "inchh")
  ic <- grep("^inchh_", names(d), value = TRUE)
  d$hh_total <- rowSums(d[ic])
  d$hh_den   <- d$hh_total - d[[INCOME_UNSPEC_COL]]
  grp <- list(sub_place = d$sp_code, main_place = d$mp_code, municipality = d$lm_code)

  shr <- function(cols) pool_share(rowSums(d[cols]), d$hh_den, grp, min_hh)
  p  <- shr(low_cols)
  d$share_low  <- p$share
  d$pool_level <- p$level
  d$share_ub   <- shr(LOW_INCOME_COLS_UB)$share
  d$share_none <- shr(NO_INCOME_COLS)$share
  d$hh_low_raw <- rowSums(d[low_cols])

  # 2. SAL polygons touching the padded grid -----------------------------------
  exdir <- file.path(tempdir(), "census2011_sal")
  shp <- file.path(exdir, "SAL_SA_2013.shp")
  if (!file.exists(shp)) {
    utils::unzip(mapdata_zip, junkpaths = TRUE, exdir = exdir,
                 files = paste0("mapdata-sa/SAL_SA_2013/SAL_SA_2013.",
                                c("shp", "shx", "dbf", "prj")))
  }
  sal <- sf::st_read(shp, quiet = TRUE, wkt_filter = sf::st_as_text(g$bb_ll)) |>
    sf::st_make_valid()
  sal$sal_code <- sprintf("%.0f", sal$SAL_CODE)
  sal <- sal[, c("sal_code", "SP_NAME", "MN_NAME")]
  names(sal)[2:3] <- c("sp_name", "lm_name")
  sal <- merge(sal, d[c("sal_code", "lm_code", "hh_total", "hh_den", "hh_low_raw",
                        "share_low", "share_ub", "share_none", "pool_level")],
               by = "sal_code", all.x = TRUE)
  if (anyNA(sal$hh_total)) stop(sum(is.na(sal$hh_total)), " SALs have no income table row.")
  sal$hh_low <- ifelse(sal$hh_den > 0, sal$share_low * sal$hh_den, 0)
  sal <- sf::st_transform(sal, terra::crs(g$tmpl))

  if (!is.null(out_gpkg)) {
    dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
    sf::st_write(sal, out_gpkg, layer = "low_income_sal2011", delete_dsn = TRUE, quiet = TRUE)
  }

  # 3. Spread on settlement, 1 km ------------------------------------------------
  settle <- settlement_weights(sanlc_file, g$sub, settlement_classes)
  v <- terra::vect(sal)
  agg <- function(f) {
    sp <- spread_households(v, g$sub, settle, f)
    list(r = terra::aggregate(sp$surface, fact = fine, fun = "sum"), n = sp$n)
  }
  low <- agg("hh_low")
  den <- agg("hh_den")

  sh   <- to_template(terra::ifel(den$r > 0, low$r / den$r, NA), g$tmpl, "low_income_share")
  outs <- write_1km(sh, out_share)
  if (!is.null(out_gpkg)) outs <- c(outs, out_gpkg)

  # QA ---------------------------------------------------------------------------
  pop <- sal$hh_den > 0
  hh_in  <- terra::global(to_template(den$r, g$tmpl, "x"), "sum", na.rm = TRUE)[[1]]
  low_in <- terra::global(to_template(low$r, g$tmpl, "x"), "sum", na.rm = TRUE)[[1]]
  vals <- terra::values(sh, mat = FALSE); vals <- vals[is.finite(vals)]
  q <- round(stats::quantile(vals, c(0.25, 0.5, 0.75)), 2)
  message(
    "SALs near the grid: ", nrow(sal), " | populated: ", sum(pop), " | pooled (<", min_hh,
    " hh): ", paste(names(table(sal$pool_level[pop])), table(sal$pool_level[pop]),
                    sep = " ", collapse = ", "),
    "\nSpread on settlement: ", den$n[["settlement"]], ", evenly: ", den$n[["even"]],
    ", by centroid: ", den$n[["centroid"]],
    "\nCRA extent: households ", round(hh_in), ", low-income ", round(low_in),
    " (", round(100 * low_in / hh_in, 1), "%)",
    "\nCells with households: ", length(vals),
    " | share p25 ", q[1], ", median ", q[2], ", p75 ", q[3]
  )
  outs
}
