#************************************************************
# 101d_download_nex_aws.R
# Fill missing NEX-GDDP-CMIP6 v2.0 files from the AWS Open Data copy
# (s3://nex-gddp-cmip6, public HTTPS) instead of the NCCS THREDDS server.
#
#   1. Inventory and availability as in 101b; invalid files are quarantined.
#   2. The v2.0 md5 index is downloaded once to the manifests folder and
#      matched to the missing year x model x variable rows.
#   3. test = TRUE fetches only n_test files and prints time per file and an
#      estimate for the full gap; test = FALSE fetches every missing file.
#      Each global file (~180-250 MB) is downloaded, md5-checked, cropped to
#      the grid of an existing valid file, written to its period folder under
#      the NCSS file name, and deleted.
#   4. Results go to manifests/aws_fetch_log.csv; 101b is then sourced to
#      rewrite the availability manifests.
#************************************************************

library(dplyr)
library(readr)
library(tidyr)
source("R/functions/nex_functions.R")

raw_dir <- "data/raw/nex_gddp_cmip6"
man_dir <- file.path(raw_dir, "manifests")
vars    <- c("pr", "tasmax", "tasmin")
use_ssps <- c("ssp126", "ssp245")
test    <- TRUE
n_test  <- 3
workers <- 2
tmp_dir <- file.path(tempdir(), "nex_aws")
dir.create(tmp_dir, showWarnings = FALSE)

models <- read_csv(file.path(man_dir, "common_complete_ensemble.csv"),
                   show_col_types = FALSE)$model
jobs <- read_csv(file.path(man_dir, "download_jobs.csv"), show_col_types = FALSE) |>
  filter(scenario %in% c("historical", use_ssps))

inv <- nex_inventory(raw_dir, jobs$period_name, vars, models)
nex_quarantine(inv, raw_dir)
inv <- inv[inv$valid, ]
av  <- nex_availability(inv, jobs, vars, models)

todo <- nex_todo(av)
if (is.null(todo)) stop("Nothing missing.")
ref <- inv$file[inv$scenario == "historical"][1]
stopifnot(length(ref) == 1, !is.na(ref))

idx <- nex_aws_index(man_dir)
message(nrow(todo), " files missing; ",
        sum(paste(todo$model, todo$scenario, todo$variable, todo$year) %in%
              paste(idx$model, idx$scenario, idx$variable, idx$year)),
        " are in the AWS v2.0 index")

if (test) todo <- todo[seq_len(min(n_test, nrow(todo))), ]

t0  <- Sys.time()
res <- nex_aws_fetch(todo, idx, raw_dir, ref, workers = workers, tmp_dir = tmp_dir)
el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

log_f <- file.path(man_dir, "aws_fetch_log.csv")
write_csv(res |> transmute(period, model, scenario, variable, year, file, mb, status,
                           time = format(Sys.time())),
          log_f, append = file.exists(log_f))
print(count(res, period, status), n = Inf)

if (test) {
  n_all <- nrow(nex_todo(av))
  ok <- res$status == "ok"
  message(sprintf("Test: %d/%d ok, %.1f min (%.2f min/file, %.0f MB/file). ",
                  sum(ok), nrow(res), el, el / nrow(res), mean(res$mb, na.rm = TRUE)),
          sprintf("Full gap of %d files: ~%.0f h and ~%.0f GB downloaded.",
                  n_all, n_all * el / nrow(res) / 60,
                  n_all * mean(res$mb, na.rm = TRUE) / 1024))
} else {
  source("R/10_download/101b_check_nex_downloads.R")
}
