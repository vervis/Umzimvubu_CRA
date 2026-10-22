#
# 101a_fix_and_resume_nex_download.R
#
# Run this BEFORE re-running the main download script.
# It corrects misplaced files from the first (partial) download run,
# then resumes downloading the missing models and SSPs.
#
# Summary of what happened in the first run:
#   - Only 4 models were downloaded (GFDL-ESM4, MPI-ESM1-2-HR,
#     MRI-ESM2-0, UKESM1-0-LL); 21 models are still missing.
#   - Only SSP2-4.5 data exists for the future periods.
#   - That ssp245 data landed in wrongly-named folders (no SSP suffix),
#     and the near-term slice ended up inside near_2021_2040_ssp126/.
#
# Steps:
#   Step 1 - Move the existing ssp245 files to their correct directories.
#   Step 2 - Re-run cmip6_dl with force = FALSE so it skips what is
#             already present and downloads everything that is missing.
#

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(stringr)
library(readr)
library(cmip6)

Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.tiff,.nc",
  VSI_CACHE = "TRUE",
  VSI_CACHE_SIZE = "500000000"
)

terraOptions(progress = 1, memfrac = 0.75)

# Region of interest (same as main script)
shp_umz <- st_read("../../Hub_boundaries/Umzimvubu_Catchment.shp", quiet = TRUE) |>
  st_make_valid() |> st_transform(4326)
shp_anzo <- st_read("../../Hub_boundaries/Alfred_Nzo.shp", quiet = TRUE) |>
  st_make_valid() |> st_transform(4326)
shp_combined <- dplyr::bind_rows(shp_umz, shp_anzo) |>
  sf::st_union() |> sf::st_make_valid() |> sf::st_as_sf()


#************************************************************
# Step 1: Fix misplaced directories from the first run
#************************************************************

relocate_nc_files <- function(from_dir, to_dir) {
  if (!dir.exists(from_dir)) {
    message("Source directory does not exist, skipping: ", from_dir)
    return(invisible(NULL))
  }
  nc_files <- list.files(from_dir, pattern = "\\.nc$", full.names = TRUE)
  if (length(nc_files) == 0) {
    message("No .nc files found in: ", from_dir)
    return(invisible(NULL))
  }
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  message("Moving ", length(nc_files), " files: ", from_dir, " -> ", to_dir)
  for (f in nc_files) {
    dest <- file.path(to_dir, basename(f))
    if (!file.exists(dest)) {
      file.rename(f, dest)
    } else {
      message("  Skipping (already exists): ", basename(f))
    }
  }
}

# a) near_2021_2040_ssp126/ holds ssp245 files — move to near_2021_2040_ssp245/
relocate_nc_files(
  from_dir = "data/raw/nex_gddp_cmip6/near_2021_2040_ssp126",
  to_dir   = "data/raw/nex_gddp_cmip6/near_2021_2040_ssp245"
)

# b) mid_2041_2070/ (no SSP suffix) holds ssp245 files
relocate_nc_files(
  from_dir = "data/raw/nex_gddp_cmip6/mid_2041_2070",
  to_dir   = "data/raw/nex_gddp_cmip6/mid_2041_2070_ssp245"
)

# c) late_2071_2100/ (no SSP suffix) holds ssp245 files
relocate_nc_files(
  from_dir = "data/raw/nex_gddp_cmip6/late_2071_2100",
  to_dir   = "data/raw/nex_gddp_cmip6/late_2071_2100_ssp245"
)

message("\nStep 1 complete. Verifying file counts...")
for (d in c("near_2021_2040_ssp245", "mid_2041_2070_ssp245", "late_2071_2100_ssp245")) {
  n <- length(list.files(file.path("data/raw/nex_gddp_cmip6", d), pattern = "\\.nc$"))
  message("  ", d, ": ", n, " .nc files")
}


#************************************************************
# Step 2: Re-run the full download (force = FALSE skips
#          files that already landed in Step 1)
#************************************************************

nex_vars      <- c("pr", "tas", "tasmin", "tasmax")
nex_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

nex_models <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CMCC-ESM2",
  "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg-LR", "FGOALS-g3",
  "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
  "KACE-1-0-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
  "MRI-ESM2-0", "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"
)

nex_jobs <- dplyr::bind_rows(
  tibble::tibble(
    period_name = "hist_1981_2014",
    scenario = "historical",
    years = list(1981:2014)
  ),
  purrr::map(
    nex_scenarios,
    function(ssp) {
      tibble::tribble(
        ~period_name,                    ~scenario, ~years,
        paste0("near_2021_2040_", ssp),  ssp,       2021:2040,
        paste0("mid_2041_2070_", ssp),   ssp,       2041:2070,
        paste0("late_2071_2100_", ssp),  ssp,       2071:2100
      )
    }
  ) |>
    purrr::list_rbind()
)

download_nex_subset <- function(period_name, scenario, years) {
  outdir <- file.path("data/raw/nex_gddp_cmip6", period_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  n_existing <- length(list.files(outdir, pattern = "\\.nc$"))
  message(
    "\nDownloading: ", period_name, " / ", scenario,
    "  (", n_existing, " files already present)"
  )

  cmip6::cmip6_dl(
    outdir    = outdir,
    aoi       = shp_combined,
    models    = nex_models,
    scenarios = scenario,
    elements  = nex_vars,
    years     = unlist(years),
    workers   = 4,
    force     = FALSE   # skips any file that already exists
  )
}

purrr::pwalk(nex_jobs, download_nex_subset)

message("\nDownload run complete.")
message("Expected file totals (full 25-model ensemble):")
message("  hist_1981_2014 : 3,400   (25 models x 4 vars x 34 years)")
message("  each near/*    : 2,000   (25 models x 4 vars x 20 years)")
message("  each mid/*     : 3,000   (25 models x 4 vars x 30 years)")
message("  each late/*    : 3,000   (25 models x 4 vars x 30 years)")
message("  GRAND TOTAL    : 35,400 files")
