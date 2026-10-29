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
# NEX-GDDP indicators are also scored per GCM: the normalisation bounds are
# fitted once to the pooled per-model layers (all GCMs x SSPs x periods) and
# applied to the ensemble median and to every GCM, so the scores are on one
# scale and individual models are not clipped by median-based bounds.
#
# Everything about an individual indicator (its file, resampling method,
# direction, bounds, weight) lives in R/30_indicators/indicator_register.csv,
# not in code. _targets.R runs make_indicator() once per active row.
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

# Per-GCM scores are stored as INT2U = round(score * GCM_SCALE); read them
# with read_gcm_scores().
GCM_SCALE <- 10000


# ----------------------------------------------------------------------------
# read_indicator_register(): load and check the register
# ----------------------------------------------------------------------------
# Returns only the active rows. `upstream` is converted to symbols so that
# tarchetypes::tar_map() makes each indicator depend on the stage-20 target
# that writes its source file. Optional `weight` (default 1) is the weight
# of the indicator within its component in stage 40.

read_indicator_register <- function(path, active_only = TRUE, syms = TRUE) {
  reg <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))

  required <- c("id", "active", "component", "subrisk", "src", "temporal",
                "resample", "direction", "normalise", "trim",
                "lower", "upper", "na_fill", "upstream")
  missing <- setdiff(required, names(reg))
  if (length(missing)) stop("Register is missing columns: ",
                            paste(missing, collapse = ", "))

  if (!"transform" %in% names(reg)) reg$transform <- NA_character_
  reg$transform <- as.character(reg$transform)
  if (!"weight" %in% names(reg)) reg$weight <- NA_real_
  reg$weight <- as.numeric(reg$weight)
  reg$weight[is.na(reg$weight)] <- 1

  if (active_only) reg <- reg[as.logical(reg$active), , drop = FALSE]

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
  bad <- reg$weight < 0
  if (any(bad)) stop("weight must be >= 0: ", paste(reg$id[bad], collapse = ", "))

  for (col in c("trim", "lower", "upper", "na_fill")) reg[[col]] <- as.numeric(reg[[col]])

  if (syms) reg$upstream <- rlang::syms(reg$upstream)
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
# gcm_sources(): per-model 201 files behind an ensemble-median `src`
# ----------------------------------------------------------------------------
# "..._{period}_ensemble_median_umz.tif" -> "..._<period>_<model>_umz.tif".
# Returns a list (one element per period) of model-named file vectors.

gcm_sources <- function(src, periods, models) {
  if (!grepl("ensemble_median", src, fixed = TRUE))
    stop("NEX src must point to an ensemble_median file: ", src)
  out <- lapply(stats::setNames(periods, periods), function(p) {
    f <- gsub("{period}", p, sub("ensemble_median", "{model}", src, fixed = TRUE), fixed = TRUE)
    stats::setNames(vapply(models, function(m) gsub("{model}", m, f, fixed = TRUE), ""), models)
  })
  missing <- unlist(out)[!file.exists(unlist(out))]
  if (length(missing)) stop("Per-model source file(s) not found:\n  ",
                            paste(utils::head(missing, 20), collapse = "\n  "), call. = FALSE)
  out
}


read_gcm_scores <- function(f) terra::rast(f) / GCM_SCALE


# ----------------------------------------------------------------------------
# make_indicator(): one register row -> value + score GeoTIFFs
# ----------------------------------------------------------------------------
# Writes (one layer per period; "baseline" for static):
#   <out_dir>/<component>/<id>_value.tif : on the template grid, real units
#   <out_dir>/<component>/<id>_score.tif : 0-1, where 1 = worst
#   <out_dir>/<component>/<id>_bounds.csv: the fitted normalisation (lo, hi
#                                          on the transformed scale)
# and for NEX indicators when `models` is given:
#   <out_dir>/<component>/gcm/<id>_score_<period>.tif : one INT2U layer per
#                                          GCM, score * GCM_SCALE
#
# Normalisation is pooled across all periods of an indicator (and, for NEX
# with `models`, across all GCMs), so a future layer can score higher than
# any baseline cell and the change stays visible.
#
#   resample  : terra::project() method. "average" when going from a finer
#               grid (fractions and 0/1 layers become proportions);
#               "bilinear" for coarse continuous fields (NEX-GDDP); "near"
#               or "mode" for classes.
#   direction : 1 if higher values = worse, -1 if higher values = better
#   normalise : "minmax" or "percentile" (rank within the pooled values)
#   trim      : optional, e.g. 0.02 clips values below the 2nd and above the
#               98th percentile before minmax, so outliers don't squash the
#               rest. Ignored when lower/upper are given.
#   lower, upper : optional fixed bounds (e.g. 0 and 1 for shares).
#   na_fill   : optional value for template cells that are NA
#               after resampling (applied before normalising).
#   transform : optional, "log1p" or "sqrt" (INDICATOR_TRANSFORMS). Applied
#               after na_fill and before normalising; trim and lower/upper
#               then refer to the transformed values. The _value.tif is
#               written untransformed.
#   models    : GCM names (201 ensemble); NULL = no per-GCM scoring.
#   upstream  : the stage-20 target that writes `src`. Not used, only there
#               so targets builds things in the right order.

make_indicator <- function(id, component, src, temporal, periods,
                           template, resample, direction,
                           normalise = "minmax", trim = NA,
                           lower = NA, upper = NA, na_fill = NA,
                           transform = NA, models = NULL,
                           out_dir = "data/processed/indicators",
                           upstream = NULL) {
  force(upstream)
  files <- indicator_sources(src, temporal, periods)
  tmpl  <- terra::rast(template)
  use_gcm <- temporal == "nex" && length(models) > 0
  has_tf <- !is.na(transform) && nzchar(transform)
  if (has_tf && !transform %in% INDICATOR_TRANSFORMS)
    stop(id, ": unknown transform '", transform, "'.")

  prep <- function(f) {
    r <- terra::rast(f)
    if (terra::nlyr(r) > 1) {
      warning(id, ": ", basename(f), " has ", terra::nlyr(r),
              " layers; using the first.", call. = FALSE)
      r <- r[[1]]
    }
    r <- terra::project(r, tmpl, method = resample)
    if (!is.na(na_fill)) r <- terra::ifel(is.na(r) & !is.na(tmpl), na_fill, r)
    terra::mask(r, tmpl)
  }
  tf <- function(r) {
    if (!has_tf) return(r)
    lo_v <- terra::global(r, "min", na.rm = TRUE)[[1]]
    if (any(is.finite(lo_v) & lo_v < 0))
      stop(id, ": transform '", transform, "' needs values >= 0 (min ",
           min(lo_v, na.rm = TRUE), ").")
    if (transform == "log1p") log1p(r) else sqrt(r)
  }

  # 1. Resample to template -------------------------------------------------
  value <- terra::rast(lapply(files, prep))
  names(value) <- names(files)
  gcm <- if (use_gcm) {
    lapply(gcm_sources(src, names(files), models), function(fs) {
      r <- terra::rast(lapply(fs, prep))
      names(r) <- names(fs)
      r
    })
  }

  # 2. Pooled values for the bounds (transformed scale) ---------------------
  pool <- if (use_gcm) unlist(lapply(gcm, function(r) terra::values(tf(r), mat = FALSE)))
          else terra::values(tf(value), mat = FALSE)
  pool <- pool[is.finite(pool)]
  if (!length(pool)) stop(id, ": no valid cells in the CRA extent.")

  # 3. Normalise ------------------------------------------------------------
  if (normalise == "percentile") {
    ref <- sort(pool)
    lo <- hi <- NA_real_
    sc <- function(r) {
      s <- terra::app(tf(r), function(x) {
        out <- rep(NA_real_, length(x))
        ok  <- is.finite(x)
        out[ok] <- findInterval(x[ok], ref) / length(ref)
        out
      })
      if (direction == -1) 1 - s else s
    }
  } else {
    lo <- if (!is.na(lower)) lower else
      if (!is.na(trim)) stats::quantile(pool, trim, names = FALSE) else min(pool)
    hi <- if (!is.na(upper)) upper else
      if (!is.na(trim)) stats::quantile(pool, 1 - trim, names = FALSE) else max(pool)
    if (hi <= lo) stop(id, ": upper bound (", hi, ") <= lower bound (", lo, ").")
    sc <- function(r) {
      s <- terra::clamp((tf(r) - lo) / (hi - lo), 0, 1, values = TRUE)
      if (direction == -1) 1 - s else s
    }
  }
  score <- sc(value)
  names(score) <- names(files)

  # 4. Write ----------------------------------------------------------------
  dir <- file.path(out_dir, component)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out <- c(
    value  = file.path(dir, paste0(id, "_value.tif")),
    score  = file.path(dir, paste0(id, "_score.tif")),
    bounds = file.path(dir, paste0(id, "_bounds.csv"))
  )
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  terra::writeRaster(value, out[["value"]], overwrite = TRUE,
                     datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(score, out[["score"]], overwrite = TRUE,
                     datatype = "FLT4S", gdal = gdal)
  utils::write.csv(data.frame(
    id = id, normalise = normalise, transform = if (has_tf) transform else NA,
    direction = direction, lo = lo, hi = hi,
    pooled = if (use_gcm) "gcm" else "layers", n_models = length(models) * use_gcm,
    n_values = length(pool)
  ), out[["bounds"]], row.names = FALSE)

  gf <- character()
  if (use_gcm) {
    gdir <- file.path(dir, "gcm")
    dir.create(gdir, showWarnings = FALSE)
    gf <- vapply(names(gcm), function(p) {
      f <- file.path(gdir, paste0(id, "_score_", p, ".tif"))
      s <- round(sc(gcm[[p]]) * GCM_SCALE)
      names(s) <- names(gcm[[p]])
      terra::writeRaster(s, f, overwrite = TRUE, datatype = "INT2U",
                         gdal = c(gdal, "PREDICTOR=2"))
      f
    }, "")
  }
  unname(c(out, gf))
}


# ----------------------------------------------------------------------------
# summarise_indicators(): QA table for all indicator files
# ----------------------------------------------------------------------------
# One row per value/score file and layer (bounds and per-GCM files are
# skipped). `coverage` is the share of template cells with a value; anything
# well below 1 means the source does not cover the whole extent (or needs
# na_fill).

summarise_indicators <- function(files, template, out) {
  files <- files[grepl("_(value|score)\\.tif$", files)]
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
