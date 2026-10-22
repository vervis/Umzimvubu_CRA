#************************************************************
# 101c_refetch_nex_failures.R
# Re-request the NEX-GDDP-CMIP6 files that 101 failed to download.
#
# The NASA NCCS THREDDS server returns intermittent HTTP 500/504 errors;
# cmip6_dl() saves these error pages as .nc files and does not retry.
#   1. Inventory and availability as in 101b; invalid files are moved to
#      data/raw/nex_gddp_cmip6/_failed/<period>/.
#   2. nex_refetch(): every missing year x model x variable is requested
#      directly from NCSS (5 immediate retries per request, then up to 30
#      rounds over the files still missing), using
#      the grid of an existing valid file; only NetCDF responses on that grid
#      are kept.
#   3. Availability is recomputed and the 101b manifests are rewritten by
#      sourcing 101b.
#************************************************************

library(dplyr)
library(readr)
library(tidyr)
source("R/functions/nex_functions.R")

raw_dir <- "data/raw/nex_gddp_cmip6"
man_dir <- file.path(raw_dir, "manifests")
vars <- c("pr", "tasmax", "tasmin")
use_ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")

models <- read_csv(file.path(man_dir, "common_complete_ensemble.csv"),
                   show_col_types = FALSE)$model
jobs <- read_csv(file.path(man_dir, "download_jobs.csv"), show_col_types = FALSE) |>
  filter(scenario %in% c("historical", use_ssps))

inv <- nex_inventory(raw_dir, jobs$period_name, vars, models)
nex_quarantine(inv, raw_dir)
inv <- inv[inv$valid, ]
av <- nex_availability(inv, jobs, vars, models)

ref <- inv$file[inv$scenario == "historical"][1]
stopifnot(length(ref) == 1, !is.na(ref))

res <- nex_refetch(av, raw_dir, ref, rounds = 30, pause = 10, workers = 4, tries = 5)
if (!is.null(res)) {
  print(res |> group_by(period, variable) |>
          summarise(requested = n(), valid = sum(valid), .groups = "drop"), n = Inf)
  write_csv(res |> filter(!valid) |> select(period, model, variable, year, url),
            file.path(man_dir, "refetch_still_missing.csv"))
}

source("R/10_download/101b_check_nex_downloads.R")
