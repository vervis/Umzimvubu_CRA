# ============================================================================
# R/functions/indicator_functions.R
#
# Stage 30 (indicators). These functions turn the stage-20 outputs (one
# variable per layer, in real units, on each source's own grid) into
# comparable CRA indicators:
#
#   1. put every layer on ONE template grid (make_template_grid(), the CRA
#      rectangle)
#   2. normalise it to 0-1 and set the direction so that 1 = worse
#
# Everything about an individual indicator (its file, resampling method,
# direction, bounds) lives in R/30_indicators/indicator_register.csv, not in
# code. _targets.R runs make_indicator() once per active row.
# ============================================================================


# ----------------------------------------------------------------------------
# make_template_grid(): the common CRA grid
# ----------------------------------------------------------------------------
# The CRA rectangle (cra_extent(), pipeline_helpers.R) at `res` m in CRA_CRS,
# every cell 1. Equal-area cells, so averages of fractions and area totals
# need no cell-area weighting.

make_template_grid <- function(catchment_shp, district_shp, out,
                               res = CRA_RES, crs = CRA_CRS) {
  template <- terra::rast(cra_extent(catchment_shp, district_shp, crs, res),
                          resolution = res, crs = crs, vals = 1L)
  names(template) <- "cra_extent"
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(template, out, overwrite = TRUE, datatype = "INT1U",
                     gdal = c("COMPRESS=DEFLATE"))
  out
}


# Transforms allowed in the register's optional `transform` column. They are
# applied to the values only for normalising (the _value.tif keeps real
# units). "log1p" = log(1 + x), for counts and densities with a long right
# tail and many zeros (zeros stay 0); "sqrt" is a milder alternative.
INDICATOR_TRANSFORMS <- c("log1p", "sqrt")


# ----------------------------------------------------------------------------
# read_indicator_register(): load and check the register
# ----------------------------------------------------------------------------
# Returns only the active rows. `upstream` is converted to symbols so that
# tarchetypes::tar_map() makes each indicator depend on the stage-20 target
# that writes its source file.

read_indicator_register <- function(path) {
  reg <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))

  required <- c("id", "active", "component", "subrisk", "src", "temporal",
                "resample", "direction", "normalise", "trim",
                "lower", "upper", "na_fill", "upstream")
  missing <- setdiff(required, names(reg))
  if (length(missing)) stop("Register is missing columns: ",
                            paste(missing, collapse = ", "))

  # Optional column: transform applied before normalising
  if (!"transform" %in% names(reg)) reg$transform <- NA_character_
  reg$transform <- as.character(reg$transform)

  reg <- reg[as.logical(reg$active), , drop = FALSE]

  if (anyDuplicated(reg$id)) stop("Duplicate indicator ids in register.")
  bad <- !is.na(reg$transform) & !reg$transform %in% INDICATOR_TRANSFORMS
  if (any(bad)) stop("transform must be empty or one of ",
                     paste(INDICATOR_TRANSFORMS, collapse = ", "), ": ",
                     paste(reg$id[bad], collapse = ", "))
  bad <- !reg$direction %in% c(1, -1)
  if (any(bad)) stop("direction must be 1 or -1: ", paste(reg$id[bad], collapse = ", "))
  bad <- !reg$temporal %in% c("static", "nex")
  if (any(bad)) stop("temporal must be static or nex: ", paste(reg$id[bad], collapse = ", "))
  bad <- !reg$normalise %in% c("minmax", "percentile")
  if (any(bad)) stop("normalise must be minmax or percentile: ",
                     paste(reg$id[bad], collapse = ", "))

  # A missing trim/bound/fill arrives as NA; keep it as a numeric NA
  for (col in c("trim", "lower", "upper", "na_fill")) reg[[col]] <- as.numeric(reg[[col]])

  reg$upstream <- rlang::syms(reg$upstream)
  tibble::as_tibble(reg)
}


# ----------------------------------------------------------------------------
# indicator_sources(): source file(s) for one indicator
# ----------------------------------------------------------------------------
# Static indicators have one file. NEX-GDDP indicators have "{period}" in
# `src`, which is replaced by each entry of `periods`. Missing files stop the
# pipeline rather than quietly dropping a scenario.

indicator_sources <- function(src, temporal, periods) {
  if (temporal == "static") {
    files <- stats::setNames(src, "baseline")
  } else {
    files <- stats::setNames(
      vapply(periods, function(p) gsub("{period}", p, src, fixed = TRUE), ""),
      periods
    )
  }
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Indicator source file(s) not found:\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
  files
}


# ----------------------------------------------------------------------------
# make_indicator(): one register row -> value + score GeoTIFFs
# ----------------------------------------------------------------------------
# Writes two multi-layer files (one layer per period; "baseline" for static):
#   <out_dir>/<component>/<id>_value.tif : on the template grid, real units
#   <out_dir>/<component>/<id>_score.tif : 0-1, where 1 = worst
#
# Normalisation is pooled across all periods of an indicator, so a future
# layer can score higher than any baseline cell and the change stays visible.
#
#   resample  : terra::project() method. "average" when going from a finer
#               grid (fractions and 0/1 layers become proportions);
#               "bilinear" for coarse continuous fields (NEX-GDDP); "near"
#               or "mode" for classes.
#   direction : 1 if higher values = worse, -1 if higher values = better
#   normalise : "minmax" or "percentile" (rank within the CRA extent)
#   trim      : optional, e.g. 0.02 clips values below the 2nd and above the
#               98th percentile before minmax, so outliers don't squash the
#               rest. Ignored when lower/upper are given.
#   lower, upper : optional fixed bounds (e.g. 0 and 1 for shares).
#               Use these where a meaningful scale exists; the score then does
#               not change if the extent changes.
#   na_fill   : optional value for template cells that are NA
#               after resampling (applied before normalising).
#   transform : optional, "log1p" or "sqrt" (INDICATOR_TRANSFORMS). Applied
#               after na_fill and before normalising; trim and lower/upper
#               then refer to the transformed values. The _value.tif is
#               written untransformed.
#   upstream  : the stage-20 target that writes `src`. Not used, only there
#               so targets builds things in the right order.

make_indicator <- function(id, component, src, temporal, periods,
                           template, resample, direction,
                           normalise = "minmax", trim = NA,
                           lower = NA, upper = NA, na_fill = NA,
                           transform = NA,
                           out_dir = "data/processed/indicators",
                           upstream = NULL) {
  force(upstream)
  files <- indicator_sources(src, temporal, periods)
  tmpl  <- terra::rast(template)

  # 1. Resample to template -------------------------------------------------
  value <- terra::rast(lapply(files, function(f) {
    r <- terra::rast(f)
    if (terra::nlyr(r) > 1) {
      warning(id, ": ", basename(f), " has ", terra::nlyr(r),
              " layers; using the first.", call. = FALSE)
      r <- r[[1]]
    }
    terra::project(r, tmpl, method = resample)
  }))
  names(value) <- names(files)

  if (!is.na(na_fill)) {
    value <- terra::ifel(is.na(value) & !is.na(tmpl), na_fill, value)
  }
  value <- terra::mask(value, tmpl)

  # 2. Optional transform (for normalising only) ----------------------------
  x <- value
  if (!is.na(transform) && nzchar(transform)) {
    if (!transform %in% INDICATOR_TRANSFORMS) {
      stop(id, ": unknown transform '", transform, "'.")
    }
    lo_v <- terra::global(value, "min", na.rm = TRUE)[[1]]
    if (any(is.finite(lo_v) & lo_v < 0)) {
      stop(id, ": transform '", transform, "' needs values >= 0 (min ",
           min(lo_v, na.rm = TRUE), ").")
    }
    x <- if (transform == "log1p") log1p(value) else sqrt(value)
  }

  # 3. Normalise (pooled across periods) -----------------------------------
  v <- terra::values(x, mat = FALSE)
  v <- v[is.finite(v)]
  if (!length(v)) stop(id, ": no valid cells in the CRA extent.")

  if (normalise == "percentile") {
    ref <- sort(v)
    score <- terra::app(x, function(x) {
      out <- rep(NA_real_, length(x))
      ok  <- is.finite(x)
      out[ok] <- findInterval(x[ok], ref) / length(ref)
      out
    })
  } else {
    lo <- if (!is.na(lower)) lower else
      if (!is.na(trim)) stats::quantile(v, trim, names = FALSE) else min(v)
    hi <- if (!is.na(upper)) upper else
      if (!is.na(trim)) stats::quantile(v, 1 - trim, names = FALSE) else max(v)
    if (hi <= lo) stop(id, ": upper bound (", hi, ") <= lower bound (", lo, ").")
    score <- terra::clamp((x - lo) / (hi - lo), 0, 1, values = TRUE)
  }

  # 4. Direction: 1 always means worse -------------------------------------
  if (direction == -1) score <- 1 - score
  names(score) <- names(files)

  # 5. Write --------------------------------------------------------------
  dir <- file.path(out_dir, component)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- c(
    value = file.path(dir, paste0(id, "_value.tif")),
    score = file.path(dir, paste0(id, "_score.tif"))
  )
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  terra::writeRaster(value, out[["value"]], overwrite = TRUE,
                     datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(score, out[["score"]], overwrite = TRUE,
                     datatype = "FLT4S", gdal = gdal)
  unname(out)
}


# ----------------------------------------------------------------------------
# summarise_indicators(): QA table for all indicator files
# ----------------------------------------------------------------------------
# One row per file and layer. `coverage` is the share of template cells
# with a value; anything well below 1 means the source does not cover the
# whole extent (or needs na_fill).

summarise_indicators <- function(files, template, out) {
  tmpl <- terra::rast(template)
  n_cells <- terra::global(tmpl, "notNA")[[1]]

  rows <- lapply(files, function(f) {
    r <- terra::rast(f)
    s <- terra::global(r, c("min", "mean", "max", "notNA"), na.rm = TRUE)
    data.frame(
      file     = basename(f),
      layer    = names(r),
      min      = s$min,
      mean     = s$mean,
      max      = s$max,
      coverage = round(s$notNA / n_cells, 3),
      row.names = NULL
    )
  })
  tab <- do.call(rbind, rows)

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(tab, out, row.names = FALSE)
  out
}
