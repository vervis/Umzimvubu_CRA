# ============================================================================
# 201_prepare_hazards_nex_gddp.R
# Prepare NEX-GDDP-CMIP6 hazard indicators for the Umzimvubu CRA
#
# Key behaviour
#   - Uses only complete model/variable/period combinations present on disk.
#   - Uses the download manifests when available; otherwise uses the same
#     25-model common ensemble and four SSPs as the download script.
#   - For SSP2-4.5, can fall back to the older unsuffixed future folders.
#   - Produces per-model indicators plus ensemble median, P10, P90,
#     inter-model SD and contributing-model count.
#   - For future periods, also produces model agreement on direction of
#     change relative to the same model's 1981-2014 value.
# ============================================================================

library(terra)
library(dplyr)
library(purrr)
library(readr)
library(lubridate)
library(SPEI)

terraOptions(progress = 1, memfrac = 0.75)

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------

raw_dir <- "data/raw/nex_gddp_cmip6"
out_root <- "data/processed/hazard"
overwrite_outputs <- TRUE
run_spei <- TRUE
cdd_threshold_mm <- 1
hot_threshold_c <- 35

nex_vars <- c("pr", "tas", "tasmin", "tasmax")
nex_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

fallback_models <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CMCC-ESM2",
  "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg-LR",
  "FGOALS-g3", "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8", "INM-CM5-0",
  "IPSL-CM6A-LR", "KACE-1-0-G", "MIROC-ES2L", "MIROC6",
  "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "NorESM2-LM",
  "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"
)

ensemble_manifest <- file.path(raw_dir, "manifests", "common_complete_ensemble.csv")
job_manifest <- file.path(raw_dir, "manifests", "download_jobs.csv")

nex_models <- if (file.exists(ensemble_manifest)) {
  read_csv(ensemble_manifest, show_col_types = FALSE)$model |> unique()
} else {
  fallback_models
}

nex_jobs <- if (file.exists(job_manifest)) {
  read_csv(job_manifest, show_col_types = FALSE) |>
    transmute(
      period_name = as.character(period_name),
      scenario = as.character(scenario),
      first_year = as.integer(first_year),
      last_year = as.integer(last_year)
    )
} else {
  bind_rows(
    tibble(
      period_name = "hist_1981_2014",
      scenario = "historical",
      first_year = 1981L,
      last_year = 2014L
    ),
    map_dfr(nex_scenarios, \(ssp) tribble(
      ~period_name,                    ~scenario, ~first_year, ~last_year,
      paste0("near_2021_2040_", ssp), ssp,       2021L,       2040L,
      paste0("mid_2041_2070_", ssp),  ssp,       2041L,       2070L,
      paste0("late_2071_2100_", ssp), ssp,       2071L,       2100L
    ))
  )
}

hist_job <- nex_jobs |> filter(scenario == "historical") |> slice(1)
if (nrow(hist_job) != 1L) stop("Exactly one historical NEX job is required.")

metric_defs <- tribble(
  ~metric,                     ~dir,                            ~prefix,
  "annual_rainfall",           "annual_rainfall",              "annual_rainfall_nex_",
  "rainfall_cv",               "rainfall_variability",         "rainfall_cv_nex_",
  "cdd",                       "consecutive_dry_days",         "cdd_mean_nex_",
  "rain_days_20mm",            "rainfall_days_20mm",           "rain_days_gt20_nex_",
  "rx1day",                    "one_day_peak_rainfall",        "rx1day_nex_",
  "peak_monthly_rainfall",     "peak_monthly_rainfall",       "peak_monthly_rainfall_nex_",
  "mean_temperature",          "mean_temperature",             "mean_temperature_nex_",
  "mean_max_temperature",      "mean_max_temperature",         "mean_max_temperature_nex_",
  "days_above_35",             "days_above_35",                "days_above_35_nex_",
  "peak_monthly_temperature",  "peak_monthly_temperature",     "peak_monthly_temperature_nex_",
  "extreme_hot_days_p99",      "extreme_hot_days_p99",         "extreme_hot_days_p99_nex_",
  "mean_min_temperature",      "mean_min_temperature",         "mean_min_temperature_nex_",
  "lowest_min_temperature",    "lowest_minimum_temperature",   "lowest_min_temperature_nex_",
  "spei12",                    "spei12",                       "spei12_mean_nex_"
)

walk(unique(metric_defs$dir), \(d)
  dir.create(file.path(out_root, d), recursive = TRUE, showWarnings = FALSE))
dir.create(file.path(out_root, "manifests"), recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------------------
# 1. Inventory downloaded NetCDF files once
# ----------------------------------------------------------------------------

all_nc <- list.files(raw_dir, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)
if (!length(all_nc)) stop("No NEX-GDDP-CMIP6 NetCDF files found under: ", raw_dir)

parse_model <- function(bn) {
  hit <- nex_models[vapply(
    nex_models,
    \(m) grepl(paste0("_day_", m, "_"), bn, fixed = TRUE),
    logical(1)
  )]
  if (length(hit)) hit[1] else NA_character_
}

raw_inventory <- tibble(
  file = all_nc,
  directory = basename(dirname(all_nc)),
  basename = basename(all_nc)
) |>
  mutate(
    variable = sub("_day_.*$", "", basename),
    model = vapply(basename, parse_model, character(1))
  )

# Parse the experiment/scenario from the filename after model name.
raw_inventory$scenario <- vapply(seq_len(nrow(raw_inventory)), function(i) {
  v <- raw_inventory$variable[i]
  m <- raw_inventory$model[i]
  b <- raw_inventory$basename[i]
  if (is.na(m)) return(NA_character_)
  prefix <- paste0(v, "_day_", m, "_")
  if (!startsWith(b, prefix)) return(NA_character_)
  sub("_.*$", "", substring(b, nchar(prefix) + 1L))
}, character(1))

raw_inventory$year <- suppressWarnings(as.integer(
  sub(".*_([0-9]{4})_v[^_]+\\.nc$", "\\1", raw_inventory$basename)
))

raw_inventory <- raw_inventory |>
  filter(variable %in% nex_vars, !is.na(model), !is.na(scenario), !is.na(year))

# Exact new folder first; old unsuffixed SSP245 folder second.
candidate_dirs <- function(job) {
  out <- job$period_name
  if (job$scenario == "ssp245") {
    out <- c(out, sub("_ssp245$", "", job$period_name))
  }
  unique(out)
}

nex_file_table <- function(job, variable, model) {
  dirs <- candidate_dirs(job)

  raw_inventory |>
    filter(
      .data$directory %in% dirs,
      .data$variable == .env$variable,
      .data$model == .env$model,
      .data$scenario == job$scenario,
      .data$year >= job$first_year,
      .data$year <= job$last_year
    ) |>
    mutate(dir_priority = match(directory, dirs)) |>
    arrange(year, dir_priority, basename) |>
    distinct(year, .keep_all = TRUE) |>
    arrange(year)
}

has_complete_nex <- function(job, variable, model) {
  x <- nex_file_table(job, variable, model)
  expected <- seq.int(job$first_year, job$last_year)
  identical(as.integer(x$year), expected)
}

available_models <- function(job, variables) {
  nex_models[vapply(nex_models, function(m) {
    all(vapply(variables, \(v) has_complete_nex(job, v, m), logical(1)))
  }, logical(1))]
}

# Availability QA.
availability_qa <- bind_rows(lapply(seq_len(nrow(nex_jobs)), function(i) {
  job <- nex_jobs[i, ]
  bind_rows(lapply(nex_vars, function(v) {
    tibble(
      period_name = job$period_name,
      scenario = job$scenario,
      variable = v,
      model = nex_models,
      n_files = vapply(nex_models, \(m) nrow(nex_file_table(job, v, m)), integer(1)),
      expected_files = job$last_year - job$first_year + 1L,
      complete = vapply(nex_models, \(m) has_complete_nex(job, v, m), logical(1))
    )
  }))
}))

write_csv(availability_qa, file.path(out_root, "manifests", "nex_processing_availability.csv"))

# ----------------------------------------------------------------------------
# 2. Raster/time helpers
# ----------------------------------------------------------------------------

read_nex_stack <- function(period_name, variable, model, scenario,
                           first_year, last_year) {
  job <- tibble(period_name, scenario, first_year, last_year)
  ft <- nex_file_table(job, variable, model)
  expected <- seq.int(first_year, last_year)

  if (!identical(as.integer(ft$year), expected)) {
    stop("Incomplete NEX files: ", period_name, " / ", variable, " / ", model,
         " (", nrow(ft), "/", length(expected), " annual files)")
  }

  r_list <- lapply(ft$file, rast)
  n_day <- vapply(r_list, \(x) as.integer(nlyr(x)), integer(1))

  calendar <- if (all(n_day == 360L)) {
    "360_day"
  } else if (all(n_day %in% c(365L, 366L))) {
    "gregorian"
  } else {
    stop("Unexpected annual layer counts for ", period_name, " / ", variable,
         " / ", model, ": ", paste(unique(n_day), collapse = ", "))
  }

  if (calendar == "gregorian") {
    expected_days <- ifelse(leap_year(ft$year), 366L, 365L)
    if (any(n_day != expected_days)) {
      stop("Gregorian layer count/year mismatch: ", period_name, " / ", variable,
           " / ", model)
    }
  }

  r <- do.call(c, r_list)
  year <- rep(ft$year, n_day)
  doy <- unlist(lapply(n_day, seq_len), use.names = FALSE)

  if (calendar == "360_day") {
    month <- ((doy - 1L) %/% 30L) + 1L
    day <- ((doy - 1L) %% 30L) + 1L
  } else {
    dates <- unlist(map2(ft$year, n_day, \(y, n)
      seq.Date(as.Date(paste0(y, "-01-01")), by = "day", length.out = n)))
    dates <- as.Date(dates, origin = "1970-01-01")
    month <- as.integer(format(dates, "%m"))
    day <- as.integer(format(dates, "%d"))
  }

  list(
    raster = r,
    index = tibble(
      layer = seq_len(nlyr(r)), year = year, day_of_year = doy,
      month = month, day_of_month = day, calendar = calendar
    ),
    files = ft$file,
    calendar = calendar
  )
}

nex_units <- function(r, variable) {
  if (variable == "pr") r <- r * 86400       # kg m-2 s-1 -> mm/day
  if (variable %in% c("tas", "tasmin", "tasmax")) r <- r - 273.15
  r
}

group_id <- function(x) match(x, unique(x))

sum_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
median_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}
min_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
max_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
sd_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}
cv_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  mu <- mean(x)
  if (!is.finite(mu) || mu == 0) NA_real_ else sd(x) / mu
}
q_na <- function(p) function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else unname(quantile(x, p, names = FALSE, type = 7))
}
count_na <- function(x) sum(is.finite(x))

annual_metric <- function(r, index, fun) {
  annual <- tapp(r, group_id(index$year), fun = fun)
  app(annual, mean_na)
}

monthly_from_daily <- function(r, index, fun = c("sum", "mean")) {
  fun <- match.arg(fun)
  ym <- paste(index$year, sprintf("%02d", index$month), sep = "-")
  groups <- group_id(ym)
  out <- tapp(r, groups, fun = if (fun == "sum") sum_na else mean_na)
  lev <- unique(ym)
  time(out) <- as.Date(paste0(lev, "-01"))
  names(out) <- paste0("x_", gsub("-", "_", lev))
  out
}

climatological_month_extreme <- function(r, index, monthly_fun, extreme_fun) {
  ym <- paste(index$year, sprintf("%02d", index$month), sep = "-")
  groups <- group_id(ym)
  monthly <- tapp(r, groups, fun = monthly_fun)
  first <- match(unique(ym), ym)
  month_no <- index$month[first]
  climatology <- tapp(monthly, group_id(month_no), fun = mean_na)
  app(climatology, extreme_fun)
}

max_cdd <- function(x, threshold = 1) {
  if (all(is.na(x))) return(NA_real_)
  dry <- !is.na(x) & x < threshold
  if (!any(dry)) return(0)
  rr <- rle(dry)
  max(rr$lengths[rr$values])
}

metric_info <- function(metric) metric_defs |> filter(.data$metric == .env$metric) |> slice(1)

metric_file <- function(metric, period_name, model) {
  z <- metric_info(metric)
  file.path(out_root, z$dir, paste0(z$prefix, period_name, "_", model, "_umz.tif"))
}

write_metric <- function(r, metric, period_name, model, layer_name) {
  names(r) <- layer_name
  writeRaster(r, metric_file(metric, period_name, model),
              overwrite = overwrite_outputs)
  invisible(r)
}

safe_run <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    warning(label, " skipped: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

# ----------------------------------------------------------------------------
# 3. Precipitation indicators
# ----------------------------------------------------------------------------

calc_pr_metrics <- function(job, model, cdd_threshold = 1) {
  message("PR: ", job$period_name, " / ", model)

  x <- read_nex_stack(job$period_name, "pr", model, job$scenario,
                      job$first_year, job$last_year)
  pr <- nex_units(x$raster, "pr")
  annual <- tapp(pr, group_id(x$index$year), fun = sum_na)

  metrics <- list(
    annual_rainfall = app(annual, mean_na),
    rainfall_cv = app(annual, cv_na),
    cdd = annual_metric(pr, x$index, \(z) max_cdd(z, cdd_threshold)),
    rain_days_20mm = annual_metric(pr > 20, x$index, sum_na),
    rx1day = annual_metric(pr, x$index, max_na),
    peak_monthly_rainfall = climatological_month_extreme(pr, x$index, sum_na, max_na)
  )

  walk2(names(metrics), metrics, \(nm, r)
    write_metric(r, nm, job$period_name, model,
                 paste0(nm, "_", job$period_name, "_", model)))

  invisible(metrics)
}

for (i in seq_len(nrow(nex_jobs))) {
  job <- nex_jobs[i, ]
  models <- available_models(job, "pr")
  message(job$period_name, " / pr: ", length(models), " complete model(s)")
  walk(models, \(m) safe_run(
    paste0("PR ", job$period_name, " / ", m),
    calc_pr_metrics(job, m, cdd_threshold_mm)
  ))
}

# ----------------------------------------------------------------------------
# 4. Temperature indicators not requiring a historical percentile threshold
# ----------------------------------------------------------------------------

calc_tas_metrics <- function(job, model) {
  message("TAS: ", job$period_name, " / ", model)
  x <- read_nex_stack(job$period_name, "tas", model, job$scenario,
                      job$first_year, job$last_year)
  r <- app(nex_units(x$raster, "tas"), mean_na)
  write_metric(r, "mean_temperature", job$period_name, model,
               paste0("mean_temperature_", job$period_name, "_", model))
}

calc_tasmax_metrics <- function(job, model, hot_threshold = 35) {
  message("TASMAX: ", job$period_name, " / ", model)
  x <- read_nex_stack(job$period_name, "tasmax", model, job$scenario,
                      job$first_year, job$last_year)
  tx <- nex_units(x$raster, "tasmax")

  metrics <- list(
    mean_max_temperature = app(tx, mean_na),
    days_above_35 = annual_metric(tx > hot_threshold, x$index, sum_na),
    peak_monthly_temperature = climatological_month_extreme(tx, x$index, mean_na, max_na)
  )

  walk2(names(metrics), metrics, \(nm, r)
    write_metric(r, nm, job$period_name, model,
                 paste0(nm, "_", job$period_name, "_", model)))
}

calc_tasmin_metrics <- function(job, model) {
  message("TASMIN: ", job$period_name, " / ", model)
  x <- read_nex_stack(job$period_name, "tasmin", model, job$scenario,
                      job$first_year, job$last_year)
  tn <- nex_units(x$raster, "tasmin")

  metrics <- list(
    mean_min_temperature = app(tn, mean_na),
    lowest_min_temperature = annual_metric(tn, x$index, min_na)
  )

  walk2(names(metrics), metrics, \(nm, r)
    write_metric(r, nm, job$period_name, model,
                 paste0(nm, "_", job$period_name, "_", model)))
}

for (i in seq_len(nrow(nex_jobs))) {
  job <- nex_jobs[i, ]

  walk(available_models(job, "tas"), \(m) safe_run(
    paste0("TAS ", job$period_name, " / ", m), calc_tas_metrics(job, m)))

  walk(available_models(job, "tasmax"), \(m) safe_run(
    paste0("TASMAX ", job$period_name, " / ", m),
    calc_tasmax_metrics(job, m, hot_threshold_c)))

  walk(available_models(job, "tasmin"), \(m) safe_run(
    paste0("TASMIN ", job$period_name, " / ", m), calc_tasmin_metrics(job, m)))
}

# ----------------------------------------------------------------------------
# 5. Extremely hot days: historical model-specific 99th percentile threshold
# ----------------------------------------------------------------------------

calc_p99_hot_days_for_model <- function(model) {
  if (!has_complete_nex(hist_job, "tasmax", model)) return(invisible(NULL))

  message("TASMAX P99 threshold: ", model)
  hx <- read_nex_stack(hist_job$period_name, "tasmax", model, hist_job$scenario,
                       hist_job$first_year, hist_job$last_year)
  htx <- nex_units(hx$raster, "tasmax")
  p99 <- app(htx, q_na(0.99))

  jobs <- nex_jobs[vapply(seq_len(nrow(nex_jobs)), \(i)
    has_complete_nex(nex_jobs[i, ], "tasmax", model), logical(1)), ]

  for (i in seq_len(nrow(jobs))) {
    job <- jobs[i, ]
    message("Extreme hot days P99: ", job$period_name, " / ", model)
    x <- read_nex_stack(job$period_name, "tasmax", model, job$scenario,
                        job$first_year, job$last_year)
    tx <- nex_units(x$raster, "tasmax")
    out <- annual_metric(tx > p99, x$index, sum_na)
    write_metric(out, "extreme_hot_days_p99", job$period_name, model,
                 paste0("extreme_hot_days_p99_", job$period_name, "_", model))
  }
}

walk(nex_models, \(m) safe_run(paste0("P99 hot days / ", m),
                                calc_p99_hot_days_for_model(m)))

# ----------------------------------------------------------------------------
# 6. SPEI-12 from precipitation - Thornthwaite PET
# ----------------------------------------------------------------------------

calc_spei_series <- function(model, jobs, write_periods) {
  read_monthly_job <- function(job, variable, fun) {
    x <- read_nex_stack(job$period_name, variable, model, job$scenario,
                        job$first_year, job$last_year)
    monthly_from_daily(nex_units(x$raster, variable), x$index, fun)
  }

  pr_list <- lapply(seq_len(nrow(jobs)), \(i) read_monthly_job(jobs[i, ], "pr", "sum"))
  tas_list <- lapply(seq_len(nrow(jobs)), \(i) read_monthly_job(jobs[i, ], "tas", "mean"))

  pr <- do.call(c, pr_list)
  tas <- do.call(c, tas_list)
  pr <- pr[[order(time(pr))]]
  tas <- tas[[order(time(tas))]]

  dates <- as.Date(time(pr))
  if (!identical(dates, as.Date(time(tas)))) stop("SPEI pr/tas dates differ for ", model)
  if (anyDuplicated(dates)) stop("Duplicated SPEI monthly dates for ", model)
  if (!is.lonlat(pr)) stop("SPEI Thornthwaite calculation requires lon/lat NEX rasters.")

  complete_dates <- seq.Date(min(dates), max(dates), by = "month")
  observed <- match(dates, complete_dates)
  ts_start <- c(year(min(complete_dates)), month(min(complete_dates)))

  pmat <- values(pr)
  tmat <- values(tas)
  lat <- xyFromCell(pr, seq_len(ncell(pr)))[, 2]
  smat <- matrix(NA_real_, nrow = ncell(pr), ncol = nlyr(pr))

  for (cell in seq_len(ncell(pr))) {
    p <- as.numeric(pmat[cell, ])
    t <- as.numeric(tmat[cell, ])
    if (all(is.na(p)) || all(is.na(t))) next

    pc <- tc <- rep(NA_real_, length(complete_dates))
    pc[observed] <- p
    tc[observed] <- t

    pet <- tryCatch(
      SPEI::thornthwaite(ts(pmax(tc, 0), start = ts_start, frequency = 12),
                         lat = lat[cell], na.rm = TRUE),
      error = \(e) rep(NA_real_, length(tc))
    )

    wb <- ts(pc - as.numeric(pet), start = ts_start, frequency = 12)
    sp <- tryCatch(
      SPEI::spei(wb, scale = 12,
                 ref.start = c(hist_job$first_year, 1),
                 ref.end = c(hist_job$last_year, 12),
                 na.rm = TRUE, verbose = FALSE)$fitted,
      error = \(e) rep(NA_real_, length(pc))
    )
    smat[cell, ] <- as.numeric(sp)[observed]
  }

  rs <- pr
  values(rs) <- smat
  time(rs) <- dates

  for (i in seq_len(nrow(write_periods))) {
    job <- write_periods[i, ]
    idx <- which(year(dates) >= job$first_year & year(dates) <= job$last_year)
    if (!length(idx)) next
    out <- app(rs[[idx]], mean_na)
    write_metric(out, "spei12", job$period_name, model,
                 paste0("spei12_", job$period_name, "_", model))
  }

  invisible(NULL)
}

if (run_spei) {
  hist_spei_models <- available_models(hist_job, c("pr", "tas"))

  # Historical SPEI once per model.
  walk(hist_spei_models, \(m) safe_run(
    paste0("SPEI historical / ", m),
    calc_spei_series(m, hist_job, hist_job)
  ))

  # Future SPEI separately for each SSP, always fitted to the same-model
  # 1981-2014 historical reference. Incomplete future periods are omitted.
  for (ssp in nex_scenarios) {
    ssp_jobs <- nex_jobs |> filter(scenario == ssp)
    for (m in hist_spei_models) {
      ok <- vapply(seq_len(nrow(ssp_jobs)), \(i)
        all(vapply(c("pr", "tas"), \(v) has_complete_nex(ssp_jobs[i, ], v, m),
                   logical(1))), logical(1))
      future_ok <- ssp_jobs[ok, ]
      if (!nrow(future_ok)) next
      jobs <- bind_rows(hist_job, future_ok)
      safe_run(paste0("SPEI ", ssp, " / ", m),
               calc_spei_series(m, jobs, future_ok))
    }
  }
}

# ----------------------------------------------------------------------------
# 7. Ensemble uncertainty and model agreement
# ----------------------------------------------------------------------------

models_for_metric <- function(job, metric) {
  switch(
    metric,
    annual_rainfall = available_models(job, "pr"),
    rainfall_cv = available_models(job, "pr"),
    cdd = available_models(job, "pr"),
    rain_days_20mm = available_models(job, "pr"),
    rx1day = available_models(job, "pr"),
    peak_monthly_rainfall = available_models(job, "pr"),
    mean_temperature = available_models(job, "tas"),
    mean_max_temperature = available_models(job, "tasmax"),
    days_above_35 = available_models(job, "tasmax"),
    peak_monthly_temperature = available_models(job, "tasmax"),
    mean_min_temperature = available_models(job, "tasmin"),
    lowest_min_temperature = available_models(job, "tasmin"),
    extreme_hot_days_p99 = intersect(
      available_models(job, "tasmax"),
      available_models(hist_job, "tasmax")
    ),
    spei12 = intersect(
      available_models(job, c("pr", "tas")),
      available_models(hist_job, c("pr", "tas"))
    ),
    character()
  )
}

write_ensemble_raster <- function(r, metric, period_name, suffix) {
  z <- metric_info(metric)
  names(r) <- paste0(metric, "_", period_name, "_", suffix)
  writeRaster(
    r,
    file.path(out_root, z$dir,
              paste0(z$prefix, period_name, "_", suffix, "_umz.tif")),
    overwrite = overwrite_outputs
  )
}

make_ensemble_products <- function(job, metric) {
  models <- models_for_metric(job, metric)
  files <- vapply(models, \(m) metric_file(metric, job$period_name, m), character(1))
  keep <- file.exists(files)
  models <- models[keep]
  files <- files[keep]

  if (!length(files)) {
    message("Ensemble skipped: ", metric, " / ", job$period_name)
    return(tibble())
  }

  message("Ensemble: ", metric, " / ", job$period_name,
          " / ", length(files), " model(s)")

  r <- rast(files)
  write_ensemble_raster(app(r, median_na), metric, job$period_name, "ensemble_median")
  write_ensemble_raster(app(r, q_na(0.10)), metric, job$period_name, "p10")
  write_ensemble_raster(app(r, q_na(0.90)), metric, job$period_name, "p90")
  write_ensemble_raster(app(r, sd_na), metric, job$period_name, "sd")
  write_ensemble_raster(app(r, count_na), metric, job$period_name, "n_models")

  paired <- character()
  if (job$scenario != "historical") {
    hist_models <- models_for_metric(hist_job, metric)
    paired <- intersect(models, hist_models)
    ff <- vapply(paired, \(m) metric_file(metric, job$period_name, m), character(1))
    hf <- vapply(paired, \(m) metric_file(metric, hist_job$period_name, m), character(1))
    ok <- file.exists(ff) & file.exists(hf)
    paired <- paired[ok]; ff <- ff[ok]; hf <- hf[ok]

    if (length(paired)) {
      change <- rast(ff) - rast(hf)

      change_direction <- app(change, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(NA_real_)
        sign(median(x))
      })

      agreement_pct <- app(change, function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(NA_real_)
        d <- sign(median(x))
        if (d == 0) 100 * mean(x == 0) else 100 * mean(sign(x) == d)
      })

      write_ensemble_raster(agreement_pct, metric, job$period_name, "agreement_pct")
      write_ensemble_raster(change_direction, metric, job$period_name, "change_direction")
      write_ensemble_raster(app(change, count_na), metric, job$period_name,
                            "agreement_n_models")
    }
  }

  tibble(
    metric = metric,
    period_name = job$period_name,
    scenario = job$scenario,
    n_models = length(models),
    models = paste(models, collapse = ";"),
    n_models_agreement = length(paired),
    agreement_models = paste(paired, collapse = ";")
  )
}

ensemble_manifest_out <- bind_rows(lapply(seq_len(nrow(nex_jobs)), function(i) {
  job <- nex_jobs[i, ]
  bind_rows(lapply(metric_defs$metric, \(metric)
    safe_run(paste0("Ensemble ", metric, " / ", job$period_name),
             make_ensemble_products(job, metric))))
}))

write_csv(ensemble_manifest_out,
          file.path(out_root, "manifests", "nex_ensemble_models_used.csv"))

message("Finished NEX-GDDP-CMIP6 hazard preparation.")
message("Availability: ", file.path(out_root, "manifests", "nex_processing_availability.csv"))
message("Ensembles:    ", file.path(out_root, "manifests", "nex_ensemble_models_used.csv"))
