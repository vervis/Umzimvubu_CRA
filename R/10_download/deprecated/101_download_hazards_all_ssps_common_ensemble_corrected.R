#
# CRA hazard-layer downloads and Umzimvubu subsetting
#
# Indicators:
#   1. Annual rainfall volume        -> NEX-GDDP-CMIP6 daily data
#   2. SPEI-12 drought severity      -> NEX-GDDP-CMIP6 daily data
#   3. Annual rainfall variability   -> NEX-GDDP-CMIP6 daily data
#   4. Consecutive dry days          -> NEX-GDDP-CMIP6 daily data
#   5. Aridity index                 -> Future Global AI/PET
#   6. Water stress                  -> WRI Aqueduct 4.0
#
# Notes:
#   - Future Global AI/PET is a ~1 km product.
#   - NEX-GDDP-CMIP6 is ~0.25 degrees, but has the daily/monthly
#     temporal structure needed for SPEI, annual variability and CDD.
#   - WRI Aqueduct is vector/sub-basin information; rasterising it to
#     the CHELSA grid does not make the native information 1 km.


#************************************************************
# 0. Packages and setup ----
#************************************************************

# install.packages(c(
#   "terra", "sf", "dplyr", "purrr", "stringr", "lubridate",
#   "readr", "curl", "SPEI", "pak", "rgee"
# ))
#
# pak::pkg_install("mt-climate-office/cmip6@*release")

library(terra)
library(sf)
library(dplyr)
library(purrr)
library(stringr)
library(lubridate)
library(readr)
library(curl)
library(SPEI)
library(cmip6)

# Make GDAL COG streaming less painful
Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.tiff,.nc",
  VSI_CACHE = "TRUE",
  VSI_CACHE_SIZE = "500000000"
)

terraOptions(progress = 1, memfrac = 0.75)

# Project paths
dir.create("data/raw/nex_gddp_cmip6", recursive = TRUE, showWarnings = FALSE)
dir.create("data/raw/ai_pet", recursive = TRUE, showWarnings = FALSE)

# Region of interest
shp_umz <- st_read("../../Hub_boundaries/Umzimvubu_Catchment.shp", quiet = TRUE) |>
  st_make_valid() |>
  st_transform(4326)
shp_anzo <- st_read("../../Hub_boundaries/Alfred_Nzo.shp", quiet = TRUE) |>
  st_make_valid() |>
  st_transform(4326)

# Combined outer extent (used for cropping rasters)
bbox_umz  <- st_bbox(shp_umz)
bbox_anzo <- st_bbox(shp_anzo)
ext_combined <- ext(
  min(bbox_umz["xmin"], bbox_anzo["xmin"]),
  max(bbox_umz["xmax"], bbox_anzo["xmax"]),
  min(bbox_umz["ymin"], bbox_anzo["ymin"]),
  max(bbox_umz["ymax"], bbox_anzo["ymax"])
)

# Create a combined sf AOI so that cmip6_dl downloads tiles
# covering BOTH the Umzimvubu catchment AND Alfred Nzo district.
shp_combined <- dplyr::bind_rows(shp_umz, shp_anzo) |>
  sf::st_union() |>
  sf::st_make_valid() |>
  sf::st_as_sf()

# Helper: crop raster to combined extent
crop_umz <- function(r) {
  r |>
    crop(ext_combined)
}



#************************************************************
# 1. NEX-GDDP-CMIP6: download spatial subsets ----
#************************************************************
# Used for:
#   - annual rainfall volume
#   - SPEI-12
#   - annual rainfall variability
#   - consecutive dry days
#
# The common complete ensemble below follows Appendix I of the
# NEX-GDDP-CMIP6 v2 Technical Note (31 May 2025). A model is retained
# only where historical plus all four Tier-1 SSPs are available for
# every variable required by this CRA:
#   pr, tas, tasmin and tasmax.
#
# Models with only partial SSP coverage, missing required variables, or
# retracted required temperature output are excluded. Availability of
# variables not used here (e.g. hurs or huss) does not affect inclusion.
#
# The cmip6 package downloads only the AOI subset via NASA NCSS/THREDDS,
# which keeps storage manageable. Historical data are downloaded once
# per model and are shared by all four future scenarios.

nex_vars <- c("pr", "tas", "tasmin", "tasmax")

nex_scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

# Common complete ensemble for historical + SSP1-2.6 + SSP2-4.5 +
# SSP3-7.0 + SSP5-8.5, across pr, tas, tasmin and tasmax.
nex_models <- c(
  "ACCESS-CM2",
  "ACCESS-ESM1-5",
  "BCC-CSM2-MR",
  "CanESM5",
  "CMCC-ESM2",
  "CNRM-CM6-1",
  "CNRM-ESM2-1",
  "EC-Earth3",
  "EC-Earth3-Veg-LR",
  "FGOALS-g3",
  "GFDL-ESM4",
  "GISS-E2-1-G",
  "INM-CM4-8",
  "INM-CM5-0",
  "IPSL-CM6A-LR",
  "KACE-1-0-G",
  "MIROC-ES2L",
  "MIROC6",
  "MPI-ESM1-2-HR",
  "MPI-ESM1-2-LR",
  "MRI-ESM2-0",
  "NorESM2-LM",
  "NorESM2-MM",
  "TaiESM1",
  "UKESM1-0-LL"
)

# Models excluded from the common complete ensemble and the relevant
# reason for this CRA's required variables.
nex_models_excluded <- tibble::tribble(
  ~model,             ~reason,
  "CESM2",            "tasmax and tasmin unavailable",
  "CESM2-WACCM",      "tasmax and tasmin unavailable; other required variables have partial SSP coverage",
  "CMCC-CM2-SR5",     "tasmax and tasmin source output retracted",
  "GFDL-CM4-gr1",     "only some SSPs available",
  "GFDL-CM4-gr2",     "only some SSPs available",
  "HadGEM3-GC31-LL",  "only some SSPs available",
  "HadGEM3-GC31-MM",  "only some SSPs available",
  "IITM-ESM",         "tasmax and tasmin unavailable; year 2100 incomplete",
  "KIOST-ESM",        "only some SSPs available",
  "NESM3",            "only some SSPs available"
)

# Download periods. Period names include the SSP so files from different
# scenarios cannot be written into the same directory.
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

# Save the exact ensemble and job manifest used by the analysis.
dir.create(
  "data/raw/nex_gddp_cmip6/manifests",
  recursive = TRUE,
  showWarnings = FALSE
)

readr::write_csv(
  tibble::tibble(
    model = nex_models,
    common_complete_ensemble = TRUE,
    required_variables = paste(nex_vars, collapse = ";"),
    scenarios = paste(nex_scenarios, collapse = ";")
  ),
  "data/raw/nex_gddp_cmip6/manifests/common_complete_ensemble.csv"
)

readr::write_csv(
  nex_models_excluded,
  "data/raw/nex_gddp_cmip6/manifests/excluded_models.csv"
)

readr::write_csv(
  nex_jobs |>
    dplyr::mutate(
      first_year = purrr::map_int(years, min),
      last_year = purrr::map_int(years, max)
    ) |>
    dplyr::select(-years),
  "data/raw/nex_gddp_cmip6/manifests/download_jobs.csv"
)

download_nex_subset <- function(period_name, scenario, years) {
  outdir <- file.path("data/raw/nex_gddp_cmip6", period_name)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message(
    "Downloading NEX-GDDP-CMIP6 subset: ",
    period_name,
    " / ", scenario,
    " / ", length(nex_models), " models",
    " / ", length(nex_vars), " variables"
  )

  cmip6::cmip6_dl(
    outdir    = outdir,
    aoi       = shp_combined,   # was shp_umz; now covers both catchments
    models    = nex_models,
    scenarios = scenario,
    elements  = nex_vars,
    years     = unlist(years),
    workers   = 4,              # safer than the default 10 for NCSS
    force     = FALSE
  )
}

purrr::pwalk(nex_jobs, download_nex_subset)


#************************************************************
# 2. Global AI/PET: aridity index ----
#************************************************************
# The Global AI/PET aridity-index dataset is distributed through ScienceDB.
# Use a URL manifest for both historical and future annual PET_AI products.
#
# Recommended manifest:
#   data/raw/ai_pet/ai_pet_aridity_manifest.csv
#
# Required columns, after standardisation:
#   period, scenario, ensemble, variable, url, tif_pattern
#
# Example rows:
#   1960-1990,historical,WorldClim_1_4,ai,<direct_zip_url>,AI.*yr.*tif
#   1970-2000,historical,WorldClim_2_1,ai,<direct_zip_url>,AI.*yr.*tif
#   2021-2040,ssp245,MultiModel_Ensemble_Consensus,ai,<direct_zip_url>,AI.*yr.*tif
#   2041-2060,ssp245,MultiModel_Ensemble_Consensus,ai,<direct_zip_url>,AI.*yr.*tif
#
# The code below downloads, unzips, crops and masks the aridity-index rasters.

standardise_ai_pet_manifest <- function(x) {

  required_cols <- c("type", "filename", "path", "url")

  if (!all(required_cols %in% names(x))) {
    stop(
      "Manifest must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }

  x |>
    dplyr::mutate(
      type     = as.character(type),
      filename = as.character(filename),
      path     = as.character(path),
      url      = as.character(url),

      # Keep only annual aridity-index products
      is_ai_annual = stringr::str_detect(path, "PET_AI_annual"),

      # Extract time period, e.g. 1970-2000, 2021-2040, 2041-2060
      period = stringr::str_extract(filename, "\\d{4}-\\d{4}"),

      # Historical rows do not have SSP codes in the filename
      .ssp_code = stringr::str_match(filename, "_(126|245|370|585)_")[, 2],
      scenario = dplyr::case_when(
        type == "historical"     ~ "historical",
        !is.na(.ssp_code)        ~ paste0("ssp", .ssp_code),
        TRUE                     ~ NA_character_        # drops unrecognised rows
      ),

      # Remove .zip
      ensemble = stringr::str_remove(filename, "\\.zip$"),

      # Future filenames:
      #   ACCESS-ESM1-5_245_2021-2040
      #   MultiModel_Ensemble_Consensus_245_2021-2040
      ensemble = stringr::str_remove(
        ensemble,
        "_(126|245|370|585)_\\d{4}-\\d{4}$"
      ),

      # Historical filenames:
      #   WorldClim_1_4_1960-1990
      #   WorldClim_2_1_1970-2000
      ensemble = stringr::str_remove(
        ensemble,
        "_\\d{4}-\\d{4}$"
      ),

      variable = "ai",

      # Keep this broad because the exact internal GeoTIFF names can vary
      tif_pattern = "AI.*\\.tif$"
    ) |>
    dplyr::filter(
      is_ai_annual,
      !is.na(period),
      !is.na(scenario),
      !is.na(ensemble),
      !is.na(url)
    ) |>
    dplyr::select(
      period,
      scenario,
      ensemble,
      variable,
      url,
      tif_pattern
    ) |>
    dplyr::distinct()
}

download_crop_ai_pet <- function(period, scenario, ensemble, variable, url, tif_pattern) {

  message("AI/PET: ", period, " / ", scenario, " / ", ensemble, " / ", variable)

  zip_file <- file.path(
    "data/raw/ai_pet",
    paste0(variable, "_", period, "_", scenario, "_", ensemble, ".zip")
  )

  unzip_dir <- file.path(
    "data/raw/ai_pet",
    paste0(variable, "_", period, "_", scenario, "_", ensemble)
  )

  out_file <- file.path(
    "data/processed/hazard/aridity_index",
    paste0(variable, "_", period, "_", scenario, "_", ensemble, "_umz.tif")
  )

  dir.create(dirname(zip_file), recursive = TRUE, showWarnings = FALSE)
  dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(zip_file)) {
    curl::curl_download(url, zip_file, mode = "wb")
  }

  zip_contents <- unzip(zip_file, list = TRUE)

  tif_inside <- zip_contents$Name[
    stringr::str_detect(
      zip_contents$Name,
      stringr::regex(tif_pattern, ignore_case = TRUE)
    )
  ]

  if (length(tif_inside) == 0) {
    stop(
      "No GeoTIFF matched tif_pattern: ", tif_pattern,
      "\nZIP file: ", zip_file
    )
  }

  if (length(tif_inside) > 1) {
    message("Multiple matching GeoTIFFs found; using first match: ", tif_inside[1])
  }

  unzip(zip_file, files = tif_inside[1], exdir = unzip_dir, overwrite = TRUE)

  tif_file <- file.path(unzip_dir, tif_inside[1])

  r <- terra::rast(tif_file) |>
    crop_umz()

  # AI values in Global AI/PET products may be integer-scaled.
  # If values look too large for an aridity index, rescale by 0.0001.
  if (tolower(variable) == "ai" && terra::global(r, "max", na.rm = TRUE)[1, 1] > 100) {
    r <- r * 0.0001
  }

  names(r) <- paste0(variable, "_", period, "_", scenario, "_", ensemble)

  terra::writeRaster(r, out_file, overwrite = TRUE)

  invisible(r)
}

manifest_file <- "data/raw/ai_pet/AI_PET_aridity_manifest.csv"

if (file.exists(manifest_file)) {

  ai_pet_manifest <- readr::read_csv(manifest_file, show_col_types = FALSE) |>
    standardise_ai_pet_manifest()

  purrr::pwalk(
    ai_pet_manifest,
    download_crop_ai_pet
  )

} else {

  message(
    "Global AI/PET aridity-index manifest not found. ",
    "Create: ", manifest_file
  )
}
