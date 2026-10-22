#************************************************************
# 101b_check_nex_downloads.R
# Check the NEX-GDDP-CMIP6 download written by 101 before running 201.
#
#   1. Inventory the v2.0 files in every period folder listed in
#      download_jobs.csv; a file is valid only if it is real NetCDF/HDF5
#      (failed requests are saved by cmip6_dl as 39-byte "FileNotFound",
#      92-byte 504, 531-byte 500 or ~9 kB THREDDS exception pages).
#   2. Write availability per period x model x variable and a per-period
#      summary to data/raw/nex_gddp_cmip6/manifests/.
#   3. If quarantine = TRUE, move invalid v2.0 files to
#      data/raw/nex_gddp_cmip6/_failed/<period>/ so that a re-run of 101
#      (force = FALSE) requests them again.
#   Only pr, tasmax and tasmin are required.
#************************************************************

library(dplyr)
library(readr)
library(tidyr)
source("R/functions/nex_functions.R")

raw_dir <- "data/raw/nex_gddp_cmip6"
man_dir <- file.path(raw_dir, "manifests")
vars <- c("pr", "tasmax", "tasmin")
quarantine <- TRUE

models <- read_csv(file.path(man_dir, "common_complete_ensemble.csv"),
                   show_col_types = FALSE)$model
jobs <- read_csv(file.path(man_dir, "download_jobs.csv"), show_col_types = FALSE)

inv <- nex_inventory(raw_dir, jobs$period_name, vars, models)
av <- nex_availability(inv, jobs, vars, models)

smry <- av |>
  group_by(period_name, scenario) |>
  summarise(
    folder_exists = dir.exists(file.path(raw_dir, first(period_name))),
    files_expected = sum(n_expected),
    files_valid = sum(n_valid),
    files_invalid = sum(n_invalid),
    files_absent = sum(n_absent),
    pct_valid = round(100 * files_valid / files_expected, 1),
    .groups = "drop"
  ) |>
  left_join(
    av |>
      group_by(period_name, model) |>
      summarise(ok = all(complete), .groups = "drop") |>
      group_by(period_name) |>
      summarise(models_complete = sum(ok),
                models_incomplete = paste(model[!ok], collapse = ";")),
    by = "period_name"
  ) |>
  arrange(match(period_name, jobs$period_name))

hist_p <- jobs$period_name[jobs$scenario == "historical"]
common <- Reduce(intersect, lapply(jobs$period_name, \(p) nex_complete(av, p)))
by_ssp <- lapply(setdiff(unique(jobs$scenario), "historical"), \(s)
  Reduce(intersect, lapply(c(hist_p, jobs$period_name[jobs$scenario == s]),
                           \(p) nex_complete(av, p))))
names(by_ssp) <- setdiff(unique(jobs$scenario), "historical")

write_csv(av, file.path(man_dir, "download_check_availability.csv"))
write_csv(smry, file.path(man_dir, "download_check_summary.csv"))
write_csv(inv |> filter(!valid) |> mutate(file = basename(file)) |> select(-valid),
          file.path(man_dir, "download_check_invalid_files.csv"))

print(smry |> select(-models_incomplete), n = Inf, width = Inf)
message("Models complete in hist + every period of each SSP: ",
        paste(names(by_ssp), lengths(by_ssp), sep = " = ", collapse = ", "))
message("Models complete in every period: ", length(common))

if (quarantine) {
  q <- nex_quarantine(inv, raw_dir)
  message("Moved ", nrow(q), " invalid files to ", file.path(raw_dir, "_failed"),
          ". Re-run 101 for the periods with files_invalid or files_absent > 0.")
}
