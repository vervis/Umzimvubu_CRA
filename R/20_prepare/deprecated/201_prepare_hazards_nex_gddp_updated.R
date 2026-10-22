#************************************************************
# 201_prepare_hazards_nex_gddp.R
# Prepare NEX-GDDP-CMIP6 hazard indicators for the CRA
#************************************************************

library(terra)
library(dplyr)
library(purrr)
library(readr)
library(SPEI)

raw_dir <- "data/raw/nex_gddp_cmip6"
out_dir <- "data/processed/hazard"
man_dir <- file.path(raw_dir, "manifests")

# Read the exact model ensemble and jobs written by 101_download_hazards...
nex_models <- read_csv(file.path(man_dir, "common_complete_ensemble.csv"),
                       show_col_types = FALSE)$model
nex_jobs <- read_csv(file.path(man_dir, "download_jobs.csv"),
                     show_col_types = FALSE) |>
  select(period_name, scenario, first_year, last_year)
nex_ssps <- setdiff(unique(nex_jobs$scenario), "historical")

out_subdirs <- c(
  "annual_rainfall", "rainfall_variability", "consecutive_dry_days",
  "rainfall_days_20mm", "one_day_peak_rainfall", "peak_monthly_rainfall",
  "extreme_hot_days", "days_above_35c", "peak_warmest_month_temperature",
  "lowest_minimum_temperature", "spei12"
)
walk(file.path(out_dir, out_subdirs), dir.create,
     recursive = TRUE, showWarnings = FALSE)


#************************************************************
# 1. Helpers ----
#************************************************************

sum_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
min_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
q99_na <- function(x) if (all(is.na(x))) NA_real_ else quantile(x, 0.99, na.rm = TRUE, names = FALSE)
cv_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || mean(x) == 0) NA_real_ else sd(x) / mean(x)
}

read_nex_stack <- function(period_name, variable, model, scenario,
                           first_year, last_year) {
  f <- list.files(file.path(raw_dir, period_name), "\\.nc$",
                  recursive = TRUE, full.names = TRUE)
  bn <- basename(f)
  keep <- grepl(paste0("(^|[_\\.-])", variable, "([_\\.-]|$)"), bn) &
    grepl(model, bn, fixed = TRUE) & grepl(scenario, bn, fixed = TRUE)
  f <- f[keep]

  if (!length(f))
    stop("No NEX files: ", period_name, " / ", variable, " / ", model)

  yr <- suppressWarnings(as.integer(sub(".*_([0-9]{4})_v.*$", "\\1", basename(f))))
  if (anyNA(yr))
    stop("Could not extract year from: ", paste(basename(f[is.na(yr)]), collapse = ", "))

  f <- f[order(yr)]; yr <- sort(yr)
  if (anyDuplicated(yr)) stop("Duplicate annual NEX files: ", model, " / ", period_name)
  expected <- seq.int(first_year, last_year)
  if (length(yr) != length(expected) || any(yr != expected))
    stop("Missing/unexpected years: ", model, " / ", period_name)

  rr <- lapply(f, rast)
  nl <- vapply(rr, nlyr, integer(1))
  cal <- if (all(nl == 360L)) "360_day" else if (all(nl %in% c(365L, 366L))) "gregorian" else "unknown"
  if (cal == "unknown")
    stop("Unexpected daily layer counts: ", paste(basename(f), nl, sep = "=", collapse = ", "))

  day <- unlist(lapply(nl, seq_len), use.names = FALSE)
  if (cal == "360_day") {
    month <- (day - 1L) %/% 30L + 1L
  } else {
    tt <- do.call(c, Map(function(r, y) {
      d <- as.Date(time(r))
      if (length(d) != nlyr(r) || anyNA(d) || any(as.integer(format(d, "%Y")) != y))
        stop("Invalid NetCDF time axis: ", model, " / ", period_name, " / ", y)
      d
    }, rr, yr))
    month <- as.integer(format(tt, "%m"))
  }

  list(
    raster = do.call(c, rr),
    index = tibble(year = rep(yr, nl), month = month, day_of_year = day),
    calendar = cal
  )
}

nex_units <- function(r, variable) {
  if (variable == "pr") r * 86400
  else if (variable %in% c("tas", "tasmin", "tasmax")) r - 273.15
  else r
}

group_id <- function(x) match(x, unique(x))

annual_metric <- function(r, idx, fun) {
  app(tapp(r, group_id(idx$year), fun = fun), mean_na)
}

monthly_from_daily <- function(r, idx, fun) {
  id <- paste(idx$year, sprintf("%02d", idx$month), sep = "-")
  x <- tapp(r, group_id(id), fun = fun)
  time(x) <- as.Date(paste0(unique(id), "-01"))
  x
}

climatological_month_extreme <- function(r, idx, monthly_fun, extreme_fun) {
  x <- monthly_from_daily(r, idx, monthly_fun)
  m <- as.integer(format(time(x), "%m"))
  app(tapp(x, group_id(m), fun = mean_na), extreme_fun)
}

max_cdd <- function(x, threshold = 1) {
  if (all(is.na(x))) return(NA_real_)
  z <- rle(!is.na(x) & x < threshold)
  if (!any(z$values)) 0 else max(z$lengths[z$values])
}

write_metric <- function(r, subdir, prefix, period_name, model, layer_name) {
  names(r) <- layer_name
  writeRaster(r,
              file.path(out_dir, subdir,
                        paste0(prefix, period_name, "_", model, "_umz.tif")),
              overwrite = TRUE)
}

ensemble_median <- function(subdir, prefix, period_name, layer_name) {
  f <- file.path(out_dir, subdir,
                 paste0(prefix, period_name, "_", nex_models, "_umz.tif"))
  if (any(!file.exists(f)))
    stop("Missing model outputs for ensemble: ", period_name, " / ", prefix)

  r <- app(rast(f), median_na)
  names(r) <- layer_name
  writeRaster(r,
              file.path(out_dir, subdir,
                        paste0(prefix, period_name, "_ensemble_median_umz.tif")),
              overwrite = TRUE)
}


#************************************************************
# 2. Precipitation hazards ----
#************************************************************
# annual_rainfall      = mean annual precipitation total (mm/year)
# rainfall_cv          = CV of annual precipitation totals
# cdd                   = mean annual maximum run with pr < 1 mm/day
# rain_days_20mm        = mean annual number of days with pr > 20 mm
# rx1day                = mean annual maximum 1-day precipitation (mm)
# peak_monthly_rainfall = wettest climatological month (mm/month)

calc_pr_metrics <- function(job, model, cdd_threshold = 1) {
  message("PR: ", job$period_name, " / ", model)
  x <- read_nex_stack(job$period_name, "pr", model, job$scenario,
                      job$first_year, job$last_year)
  pr <- nex_units(x$raster, "pr")
  yr <- group_id(x$index$year)
  annual <- tapp(pr, yr, fun = sum_na)

  metrics <- list(
    annual_rainfall = app(annual, mean_na),
    rainfall_cv = app(annual, cv_na),
    cdd = annual_metric(pr, x$index, function(z) max_cdd(z, cdd_threshold)),
    rain_days_20mm = annual_metric(pr > 20, x$index, sum_na),
    rx1day = annual_metric(pr, x$index, max_na),
    peak_monthly_rainfall = climatological_month_extreme(pr, x$index, sum_na, max_na)
  )

  write_metric(metrics$annual_rainfall, "annual_rainfall", "annual_rainfall_nex_",
               job$period_name, model, paste0("annual_rainfall_", job$period_name, "_", model))
  write_metric(metrics$rainfall_cv, "rainfall_variability", "rainfall_cv_nex_",
               job$period_name, model, paste0("rainfall_cv_", job$period_name, "_", model))
  write_metric(metrics$cdd, "consecutive_dry_days", "cdd_mean_nex_",
               job$period_name, model, paste0("cdd_mean_", job$period_name, "_", model))
  write_metric(metrics$rain_days_20mm, "rainfall_days_20mm", "rain_days_gt20_nex_",
               job$period_name, model, paste0("rain_days_gt20_", job$period_name, "_", model))
  write_metric(metrics$rx1day, "one_day_peak_rainfall", "rx1day_nex_",
               job$period_name, model, paste0("rx1day_", job$period_name, "_", model))
  write_metric(metrics$peak_monthly_rainfall, "peak_monthly_rainfall", "peak_monthly_rainfall_nex_",
               job$period_name, model, paste0("peak_monthly_rainfall_", job$period_name, "_", model))
}

for (i in seq_len(nrow(nex_jobs)))
  for (m in nex_models)
    calc_pr_metrics(nex_jobs[i, ], m)


#************************************************************
# 3. Temperature hazards ----
#************************************************************
# extreme_hot_days = mean annual days above each model/cell's historical
#                    (1981-2014) 99th-percentile daily Tmax threshold
# days_above_35c    = mean annual days with Tmax > 35 C
# warmest_month     = highest climatological monthly mean Tmax
# lowest_tmin       = mean annual minimum daily Tmin

calc_temperature_metrics <- function(job, model, tx99, tmax_data = NULL) {
  message("Temperature: ", job$period_name, " / ", model)

  xmax <- if (is.null(tmax_data))
    read_nex_stack(job$period_name, "tasmax", model, job$scenario,
                   job$first_year, job$last_year)
  else tmax_data

  xmin <- read_nex_stack(job$period_name, "tasmin", model, job$scenario,
                         job$first_year, job$last_year)
  if (!identical(xmax$index, xmin$index))
    stop("tasmax/tasmin dates do not align: ", model, " / ", job$period_name)

  tmax <- nex_units(xmax$raster, "tasmax")
  tmin <- nex_units(xmin$raster, "tasmin")

  metrics <- list(
    extreme_hot_days = annual_metric(tmax > tx99, xmax$index, sum_na),
    days_above_35c = annual_metric(tmax > 35, xmax$index, sum_na),
    warmest_month = climatological_month_extreme(tmax, xmax$index, mean_na, max_na),
    lowest_tmin = annual_metric(tmin, xmin$index, min_na)
  )

  write_metric(metrics$extreme_hot_days, "extreme_hot_days", "extreme_hot_days_nex_",
               job$period_name, model, paste0("extreme_hot_days_", job$period_name, "_", model))
  write_metric(metrics$days_above_35c, "days_above_35c", "days_gt35c_nex_",
               job$period_name, model, paste0("days_gt35c_", job$period_name, "_", model))
  write_metric(metrics$warmest_month, "peak_warmest_month_temperature", "warmest_month_tmax_nex_",
               job$period_name, model, paste0("warmest_month_tmax_", job$period_name, "_", model))
  write_metric(metrics$lowest_tmin, "lowest_minimum_temperature", "lowest_tmin_nex_",
               job$period_name, model, paste0("lowest_tmin_", job$period_name, "_", model))
}

hist_job <- nex_jobs |> filter(scenario == "historical")
if (nrow(hist_job) != 1) stop("Expected exactly one historical NEX job")

for (m in nex_models) {
  hist_tmax <- read_nex_stack(hist_job$period_name, "tasmax", m, "historical",
                              hist_job$first_year, hist_job$last_year)
  tx99 <- app(nex_units(hist_tmax$raster, "tasmax"), q99_na)

  for (i in seq_len(nrow(nex_jobs))) {
    j <- nex_jobs[i, ]
    calc_temperature_metrics(
      j, m, tx99,
      tmax_data = if (j$scenario == "historical") hist_tmax else NULL
    )
  }
}


#************************************************************
# 4. SPEI-12 ----
#************************************************************
# Monthly precipitation - Thornthwaite PET, standardised against 1981-2014.
# Each SSP is calculated separately using the same model-specific baseline.

read_monthly <- function(job, variable, model, fun) {
  x <- read_nex_stack(job$period_name, variable, model, job$scenario,
                      job$first_year, job$last_year)
  monthly_from_daily(nex_units(x$raster, variable), x$index, fun)
}

calc_spei12 <- function(model, ssp, pr_hist, tas_hist, write_historical = FALSE) {
  message("SPEI-12: ", ssp, " / ", model)

  jobs <- nex_jobs |> filter(scenario == ssp)
  pr <- c(pr_hist, do.call(c, map(seq_len(nrow(jobs)),
                                ~ read_monthly(jobs[.x, ], "pr", model, sum_na))))
  tas <- c(tas_hist, do.call(c, map(seq_len(nrow(jobs)),
                                  ~ read_monthly(jobs[.x, ], "tas", model, mean_na))))

  pr <- pr[[order(time(pr))]]
  tas <- tas[[order(time(tas))]]
  dates <- as.Date(time(pr))
  if (anyNA(dates) || anyDuplicated(dates) || !identical(dates, as.Date(time(tas))))
    stop("Monthly pr/tas dates do not align: ", model, " / ", ssp)

  # 2015-2020 were not downloaded. Insert them as NA so the monthly time
  # series remains correctly indexed and SPEI-12 does not bridge the gap.
  all_dates <- seq.Date(min(dates), max(dates), by = "month")
  obs <- match(dates, all_dates)
  pr_mat <- values(pr)
  tas_mat <- values(tas)

  if (!is.lonlat(pr)) stop("SPEI Thornthwaite requires lon/lat NEX rasters")
  lat <- xyFromCell(pr, seq_len(ncell(pr)))[, 2]
  ans <- matrix(NA_real_, ncell(pr), nlyr(pr))
  ts_start <- c(as.integer(format(min(all_dates), "%Y")),
                as.integer(format(min(all_dates), "%m")))

  for (cell in seq_len(ncell(pr))) {
    p <- t <- rep(NA_real_, length(all_dates))
    p[obs] <- pr_mat[cell, ]
    t[obs] <- tas_mat[cell, ]
    if (all(is.na(p)) || all(is.na(t))) next

    pet <- tryCatch(
      as.numeric(SPEI::thornthwaite(
        ts(t, start = ts_start, frequency = 12),
        lat = lat[cell], na.rm = TRUE, verbose = FALSE
      )),
      error = function(e) rep(NA_real_, length(t))
    )

    sp <- tryCatch(
      SPEI::spei(
        ts(p - pet, start = ts_start, frequency = 12),
        scale = 12, ref.start = c(1981, 1), ref.end = c(2014, 12),
        na.rm = TRUE, verbose = FALSE
      )$fitted,
      error = function(e) rep(NA_real_, length(t))
    )
    ans[cell, ] <- as.numeric(sp)[obs]
  }

  spei <- pr
  values(spei) <- ans
  names(spei) <- paste0("spei12_", format(dates, "%Y_%m"))
  time(spei) <- dates
  writeRaster(spei,
              file.path(out_dir, "spei12",
                        paste0("spei12_monthly_nex_", ssp, "_", model, "_umz.tif")),
              overwrite = TRUE)

  jobs_out <- if (write_historical) bind_rows(hist_job, jobs) else jobs
  for (i in seq_len(nrow(jobs_out))) {
    j <- jobs_out[i, ]
    use <- which(dates >= as.Date(paste0(j$first_year, "-01-01")) &
                 dates <= as.Date(paste0(j$last_year, "-12-31")))
    write_metric(app(spei[[use]], mean_na), "spei12", "spei12_mean_nex_",
                 j$period_name, model,
                 paste0("spei12_mean_", j$period_name, "_", model))
  }
}

for (m in nex_models) {
  pr_hist <- read_monthly(hist_job, "pr", m, sum_na)
  tas_hist <- read_monthly(hist_job, "tas", m, mean_na)
  for (s in nex_ssps)
    calc_spei12(m, s, pr_hist, tas_hist,
                write_historical = identical(s, nex_ssps[1]))
}


#************************************************************
# 5. Ensemble medians ----
#************************************************************

ensemble_specs <- tribble(
  ~subdir,                           ~prefix,                       ~layer_prefix,
  "annual_rainfall",                "annual_rainfall_nex_",        "annual_rainfall_",
  "rainfall_variability",           "rainfall_cv_nex_",            "rainfall_cv_",
  "consecutive_dry_days",           "cdd_mean_nex_",               "cdd_mean_",
  "rainfall_days_20mm",             "rain_days_gt20_nex_",         "rain_days_gt20_",
  "one_day_peak_rainfall",          "rx1day_nex_",                 "rx1day_",
  "peak_monthly_rainfall",          "peak_monthly_rainfall_nex_",  "peak_monthly_rainfall_",
  "extreme_hot_days",               "extreme_hot_days_nex_",       "extreme_hot_days_",
  "days_above_35c",                 "days_gt35c_nex_",             "days_gt35c_",
  "peak_warmest_month_temperature", "warmest_month_tmax_nex_",     "warmest_month_tmax_",
  "lowest_minimum_temperature",     "lowest_tmin_nex_",            "lowest_tmin_",
  "spei12",                         "spei12_mean_nex_",             "spei12_mean_"
)

for (p in nex_jobs$period_name)
  for (i in seq_len(nrow(ensemble_specs)))
    ensemble_median(
      ensemble_specs$subdir[i], ensemble_specs$prefix[i], p,
      paste0(ensemble_specs$layer_prefix[i], p, "_ensemble_median")
    )

#************************************************************
# End ----
#************************************************************
