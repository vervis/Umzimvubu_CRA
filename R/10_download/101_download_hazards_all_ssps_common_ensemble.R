# ============================================================================
# 101_download_hazards_all_ssps_common_ensemble.R
#
# Downloads daily NEX-GDDP-CMIP6 v2 subsets (pr, tasmin, tasmax) for the CRA
# rectangle (cra_roi(): WGS84 bbox of the 1 km template, i.e. of the
# Umzimvubu catchment + Alfred Nzo) via NASA NCSS/THREDDS (cmip6 package).
#
# Structure
#   1. Setup and ROI.
#   2. Common complete ensemble: 25 models with historical + SSP1-2.6,
#      SSP2-4.5, SSP3-7.0 and SSP5-8.5 for every required variable
#      (NEX-GDDP-CMIP6 v2 Technical Note, Appendix I), plus excluded models.
#   3. Jobs: hist 1981-2014 and near/mid/late periods per SSP. Manifests are
#      written to data/raw/nex_gddp_cmip6/manifests/.
#   4. Download. Existing files are skipped (force = FALSE); run 101b to
#      check the downloads and move failed files out before re-running.
# ============================================================================

library(sf)
library(dplyr)
library(purrr)
library(readr)
library(cmip6)

source("R/functions/pipeline_helpers.R")

Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.tiff,.nc",
  VSI_CACHE = "TRUE",
  VSI_CACHE_SIZE = "500000000"
)

raw_dir <- "data/raw/nex_gddp_cmip6"
man_dir <- file.path(raw_dir, "manifests")
dir.create(man_dir, recursive = TRUE, showWarnings = FALSE)

# cmip6_dl needs an sf object with at least one attribute column
aoi <- st_sf(id = 1L, geometry = st_as_sfc(cra_roi()))


# 2. Ensemble ------------------------------------------------------------------

nex_vars <- c("pr", "tasmin", "tasmax")
nex_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

nex_models <- c(
  "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CanESM5", "CMCC-ESM2",
  "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "EC-Earth3-Veg-LR", "FGOALS-g3",
  "GFDL-ESM4", "GISS-E2-1-G", "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR",
  "KACE-1-0-G", "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR",
  "MRI-ESM2-0", "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"
)

nex_excluded <- tibble::tribble(
  ~model,            ~reason,
  "CESM2",           "tasmax and tasmin unavailable",
  "CESM2-WACCM",     "tasmax and tasmin unavailable; partial SSP coverage",
  "CMCC-CM2-SR5",    "tasmax and tasmin source output retracted",
  "GFDL-CM4-gr1",    "only some SSPs available",
  "GFDL-CM4-gr2",    "only some SSPs available",
  "HadGEM3-GC31-LL", "only some SSPs available",
  "HadGEM3-GC31-MM", "only some SSPs available",
  "IITM-ESM",        "tasmax and tasmin unavailable; year 2100 incomplete",
  "KIOST-ESM",       "only some SSPs available",
  "NESM3",           "only some SSPs available"
)


# 3. Jobs and manifests --------------------------------------------------------

nex_jobs <- bind_rows(
  tibble(period_name = "hist_1981_2014", scenario = "historical", years = list(1981:2014)),
  map(nex_scenarios, function(s) tibble(
    period_name = paste0(c("near_2021_2040_", "mid_2041_2070_", "late_2071_2100_"), s),
    scenario    = s,
    years       = list(2021:2040, 2041:2070, 2071:2100)
  )) |> bind_rows()
)

write_csv(
  tibble(model = nex_models, common_complete_ensemble = TRUE,
         required_variables = paste(nex_vars, collapse = ";"),
         scenarios = paste(nex_scenarios, collapse = ";")),
  file.path(man_dir, "common_complete_ensemble.csv")
)
write_csv(nex_excluded, file.path(man_dir, "excluded_models.csv"))
write_csv(
  nex_jobs |>
    mutate(first_year = map_int(years, min), last_year = map_int(years, max)) |>
    select(-years),
  file.path(man_dir, "download_jobs.csv")
)


# 4. Download ------------------------------------------------------------------

download_nex_subset <- function(period_name, scenario, years) {
  od <- file.path(raw_dir, period_name)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  message("NEX-GDDP-CMIP6: ", period_name, " / ", length(nex_models), " models")
  cmip6_dl(
    outdir    = od,
    aoi       = aoi,
    models    = nex_models,
    scenarios = scenario,
    elements  = nex_vars,
    years     = unlist(years),
    latest    = FALSE,
    workers   = 4,
    force     = FALSE
  )
}

pwalk(nex_jobs, download_nex_subset)
chk <- readr::read_csv("data/raw/nex_gddp_cmip6/manifests/download_check_summary.csv",
                       show_col_types = FALSE)
redo <- chk$period_name[chk$files_invalid + chk$files_absent > 0]
purrr::pwalk(dplyr::filter(nex_jobs, period_name %in% redo), download_nex_subset)
