#************************************************************
# 201_prepare_hazards_nex_gddp.R
# NEX-GDDP-CMIP6 hazard indicators for SR1 (water insecurity)
#
# Indicators (on the NEX grid of the 101 download, i.e. the CRA rectangle):
#   annual_rainfall  mean annual precipitation (mm/yr)
#   rainfall_cv      CV of annual precipitation totals
#   cdd              mean annual maximum run of days with pr < 1 mm
#   aridity_index    mean annual P / mean annual PET, PET by Hargreaves
#                    (FAO-56 eq. 52, daily, from tasmax/tasmin + latitude);
#                    annual PET is written as a component layer.
#
#   0. Settings, ensemble manifest and jobs from 101.
#   1. Inventory: one valid v2.0 NetCDF per variable/model/year
#      (error pages saved as .nc are ignored); availability manifest.
#   2. Ensemble membership: ens_scope "common" = models complete for hist
#      and every processed period (one consistent ensemble, as stage 30
#      normalises across periods); "period" = complete for hist + period.
#   3. Per-model metrics for every model complete in a period.
#   4. Ensemble median, p10, p90, sd, n_models; for future periods the
#      share of models agreeing on the sign of change vs their own hist.
#   5. Manifests: nex_input_availability.csv and nex_ensemble_models_used.csv
#      (the latter is the file _targets.R tracks).
#************************************************************

library(terra)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
source("R/functions/nex_functions.R")

raw_dir <- "data/raw/nex_gddp_cmip6"
out_dir <- "data/processed/hazard"
man_dir <- file.path(out_dir, "manifests")
vars <- c("pr", "tasmax", "tasmin")
cdd_thr <- 1
ens_scope <- "common"
min_models <- 10
use_ssps <- c("ssp126", "ssp245")   # same as NEX_SSPS in _targets.R
use_periods <- NULL   # NULL = hist + every period of use_ssps whose folder exists
redo <- TRUE

mets <- tribble(
  ~id,               ~dir,                   ~prefix,
  "annual_rainfall", "annual_rainfall",      "annual_rainfall_nex_",
  "rainfall_cv",     "rainfall_variability", "rainfall_cv_nex_",
  "cdd",             "consecutive_dry_days", "cdd_mean_nex_",
  "pet",             "pet_hargreaves",       "pet_hargreaves_nex_",
  "aridity_index",   "aridity_index",        "aridity_index_nex_"
)
walk(file.path(out_dir, c(mets$dir, "manifests")), dir.create,
     recursive = TRUE, showWarnings = FALSE)

models <- read_csv(file.path(raw_dir, "manifests", "common_complete_ensemble.csv"),
                   show_col_types = FALSE)$model
jobs <- read_csv(file.path(raw_dir, "manifests", "download_jobs.csv"),
                 show_col_types = FALSE) |>
  select(period_name, scenario, first_year, last_year)
hist_p <- jobs$period_name[jobs$scenario == "historical"]
if (length(hist_p) != 1) stop("Expected exactly one historical job")

if (is.null(use_periods))
  use_periods <- jobs$period_name[jobs$scenario %in% use_ssps &
                                   dir.exists(file.path(raw_dir, jobs$period_name))]
use_periods <- union(hist_p, use_periods)
jobs <- jobs[jobs$period_name %in% use_periods, ]


#************************************************************
# 1. Inventory ----
#************************************************************

inv <- nex_inventory(raw_dir, jobs$period_name, vars, models)
av <- nex_availability(inv, jobs, vars, models)
write_csv(av, file.path(man_dir, "nex_input_availability.csv"))

done <- setNames(lapply(jobs$period_name, \(p) nex_complete(av, p)), jobs$period_name)
message("Complete models per period:\n",
        paste0("  ", names(done), ": ", lengths(done), collapse = "\n"))


#************************************************************
# 2. Ensemble membership ----
#************************************************************

ens <- if (ens_scope == "common") {
  cm <- Reduce(intersect, done)
  setNames(rep(list(cm), length(done)), names(done))
} else {
  lapply(setNames(names(done), names(done)), \(p)
    if (p == hist_p) done[[p]] else intersect(done[[p]], done[[hist_p]]))
}

n_ens <- lengths(ens)
if (any(n_ens < min_models))
  stop("Fewer than ", min_models, " models in the ensemble for: ",
       paste(names(n_ens)[n_ens < min_models], collapse = ", "),
       ". See ", file.path(man_dir, "nex_input_availability.csv"),
       "; finish the download (101, 101b), drop periods via use_periods, ",
       "or set ens_scope = \"period\".")


#************************************************************
# 3. Per-model metrics ----
#************************************************************

out_file <- function(i, period, model)
  file.path(out_dir, mets$dir[i], paste0(mets$prefix[i], period, "_", model, "_umz.tif"))

read_var <- function(period, model, var, yrs) {
  x <- inv[inv$valid & inv$period == period & inv$model == model &
             inv$variable == var, ]
  x <- x[match(yrs, x$year), ]
  rr <- lapply(x$file, rast)
  nl <- as.integer(vapply(rr, nlyr, numeric(1)))
  if (!all(nl %in% c(360L, 365L, 366L)))
    stop("Unexpected daily layer count: ", var, " / ", model, " / ", period)
  for (r in rr[-1]) if (!compareGeom(rr[[1]], r, stopOnError = FALSE))
    stop("Grids differ within ", var, " / ", model, " / ", period)
  doy <- unlist(lapply(nl, \(n) if (n == 360L) (seq_len(n) - 0.5) * 365.25 / 360 + 0.5
                       else seq_len(n)))
  list(grid = rr[[1]][[1]], v = do.call(cbind, lapply(rr, values)),
       year = rep(yrs, nl), doy = doy)
}

ra_fao56 <- function(lat, doy) {
  phi <- lat * pi / 180
  dr <- 1 + 0.033 * cos(2 * pi * doy / 365)
  dec <- 0.409 * sin(2 * pi * doy / 365 - 1.39)
  ws <- acos(pmin(pmax(-tan(phi) %o% tan(dec), -1), 1))
  24 * 60 / pi * 0.082 * rep(dr, each = length(phi)) *
    (ws * (sin(phi) %o% sin(dec)) + (cos(phi) %o% cos(dec)) * sin(ws))
}

annual <- function(m, yr) t(rowsum(t(m), yr, reorder = FALSE))

max_run <- function(x) {
  if (anyNA(x)) return(NA_real_)
  z <- rle(x)
  if (any(z$values)) max(z$lengths[z$values]) else 0
}

calc_model <- function(period, model) {
  if (!redo && all(file.exists(map_chr(seq_len(nrow(mets)), out_file, period, model))))
    return(invisible())
  message("201: ", period, " / ", model)
  j <- jobs[jobs$period_name == period, ]
  yrs <- seq.int(j$first_year, j$last_year)

  p <- read_var(period, model, "pr", yrs)
  tx <- read_var(period, model, "tasmax", yrs)
  tn <- read_var(period, model, "tasmin", yrs)
  if (!identical(p$year, tx$year) || !identical(p$year, tn$year))
    stop("pr/tasmax/tasmin days do not align: ", model, " / ", period)
  if (!compareGeom(p$grid, tx$grid, stopOnError = FALSE) ||
      !compareGeom(p$grid, tn$grid, stopOnError = FALSE))
    stop("pr/tasmax/tasmin grids differ: ", model, " / ", period)
  if (!is.lonlat(p$grid)) stop("NEX grid is not lon/lat: ", model)

  pr <- p$v * 86400
  tmax <- tx$v - 273.15
  tmin <- tn$v - 273.15
  lat <- yFromCell(p$grid, seq_len(ncell(p$grid)))

  pa <- annual(pr, p$year)
  pet <- 0.0023 * 0.408 * ra_fao56(lat, p$doy) *
    ((tmax + tmin) / 2 + 17.8) * sqrt(pmax(tmax - tmin, 0))
  pet_a <- annual(pet, p$year)
  dry <- pr < cdd_thr
  cdd <- sapply(yrs, \(y) apply(dry[, p$year == y, drop = FALSE], 1, max_run))

  map_p <- rowMeans(pa)
  map_pet <- rowMeans(pet_a)
  vals <- list(
    annual_rainfall = map_p,
    rainfall_cv = ifelse(map_p > 0, apply(pa, 1, sd) / map_p, NA_real_),
    cdd = rowMeans(cdd),
    pet = map_pet,
    aridity_index = ifelse(map_pet > 0, map_p / map_pet, NA_real_)
  )

  for (i in seq_len(nrow(mets))) {
    r <- rast(p$grid)
    values(r) <- vals[[mets$id[i]]]
    names(r) <- paste0(sub("_nex_$", "_", mets$prefix[i]), period, "_", model)
    writeRaster(r, out_file(i, period, model), overwrite = TRUE)
  }
}

for (pp in jobs$period_name) for (m in done[[pp]]) calc_model(pp, m)


#************************************************************
# 4. Ensemble statistics and model agreement ----
#************************************************************

write_layer <- function(r, i, period, stat) {
  names(r) <- paste0(sub("_nex_$", "_", mets$prefix[i]), period, "_", stat)
  writeRaster(r, file.path(out_dir, mets$dir[i],
                           paste0(mets$prefix[i], period, "_", stat, "_umz.tif")),
              overwrite = TRUE)
}

ens_summary <- function(i, period) {
  ms <- ens[[period]]
  r <- rast(map_chr(ms, \(m) out_file(i, period, m)))
  st <- app(r, \(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(c(NA, NA, NA, NA, 0))
    c(median(x), quantile(x, c(0.1, 0.9), names = FALSE),
      if (length(x) > 1) sd(x) else NA, length(x))
  })
  stn <- c("median", "p10", "p90", "sd", "n_models")
  for (k in seq_along(stn)) write_layer(st[[k]], i, period, paste0("ensemble_", stn[k]))
  if (period == hist_p) return(invisible())

  h <- rast(map_chr(ms, \(m) out_file(i, hist_p, m)))
  n <- length(ms)
  ag <- app(c(r, h), \(x) {
    d <- x[seq_len(n)] - x[n + seq_len(n)]
    d <- d[is.finite(d)]
    if (!length(d)) return(c(NA, NA))
    s <- sign(median(d))
    c(100 * mean(sign(d) == s), s)
  })
  write_layer(ag[[1]], i, period, "agreement_pct")
  write_layer(ag[[2]], i, period, "change_direction")
}

for (pp in jobs$period_name) for (i in seq_len(nrow(mets))) ens_summary(i, pp)


#************************************************************
# 5. Manifest ----
#************************************************************

ens_out <- imap_dfr(ens, \(ms, pp) tibble(
  period_name = pp,
  scenario = jobs$scenario[jobs$period_name == pp],
  model = ms,
  ens_scope = ens_scope
))
write_csv(ens_out, file.path(man_dir, "nex_ensemble_models_used.csv"))
message("201 done: ", nrow(jobs), " periods, ensemble sizes ",
        paste(range(n_ens), collapse = "-"))
