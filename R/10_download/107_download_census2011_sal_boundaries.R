# ============================================================================
# 107_download_census2011_sal_boundaries.R
#
# Census 2011 Small Area Layer (SAL) boundaries, from DataFirst (UCT):
#   South African Census Community Profiles 2011 (zaf-statssa-saccp-2011-v1)
#   https://www.datafirst.uct.ac.za/dataportal/index.php/catalog/517
#
# MANUAL STEP: DataFirst requires a (free) login, so the download cannot be
# scripted. Once only:
#   1. Log in at https://www.datafirst.uct.ac.za/dataportal/index.php/auth/login
#   2. Go to https://www.datafirst.uct.ac.za/dataportal/index.php/catalog/517/get-microdata
#      and accept the terms.
#   3. Download the spatial / GIS (Small Area Layer) file.
#   4. Save it, unchanged, in data/raw/census2011/datafirst_517/
#      (a .zip, or an already-extracted folder containing a .shp or .gdb).
# Then run this script. It finds the SAL layer, keeps the SALs that
# intersect the CRA ROI, and writes a clean GeoPackage.
#
# Outputs (data/raw/boundaries/census2011_sal/):
#   census2011_sal_cra.gpkg        SALs intersecting the ROI, EPSG:4326
#   census2011_sal_manifest.csv    source file, checksum, feature count
#   census2011_sal.url             DataFirst catalogue link
#
# Standardised fields: sal_code (character), plus the original fields.
#
# _targets.R entry (10 Download section):
#   download_target(
#     "raw_sal_2011",
#     "R/10_download/107_download_census2011_sal_boundaries.R",
#     "data/raw/boundaries/census2011_sal/census2011_sal_cra.gpkg"
#   ),
# ============================================================================

library(sf)

# ============================================================================
# 1. Paths ####
# ============================================================================

src_dir <- "data/raw/census2011/datafirst_517"
out_dir <- "data/raw/boundaries/census2011_sal"
out     <- file.path(out_dir, "census2011_sal_cra.gpkg")

catchment_shp <- "data/raw/boundaries/Umzimvubu_Catchment.shp"
district_shp  <- "data/raw/boundaries/Alfred_Nzo.shp"

catalogue_url <- "https://www.datafirst.uct.ac.za/dataportal/index.php/catalog/517"

dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================================
# 2. Locate the downloaded SAL layer ####
# ============================================================================

zips <- list.files(src_dir, pattern = "\\.zip$", full.names = TRUE,
                   ignore.case = TRUE)

# Extract any zips into src_dir/extracted/<zip name>/ (only once)
for (z in zips) {
  exdir <- file.path(src_dir, "extracted",
                     tools::file_path_sans_ext(basename(z)))
  if (!dir.exists(exdir)) {
    message("Extracting ", basename(z))
    unzip(z, exdir = exdir)
  }
}

shps <- list.files(src_dir, pattern = "\\.shp$", recursive = TRUE,
                   full.names = TRUE, ignore.case = TRUE)
gdbs <- list.dirs(src_dir, recursive = TRUE)
gdbs <- gdbs[grepl("\\.gdb$", gdbs, ignore.case = TRUE)]

candidates <- c(shps, gdbs)
if (length(candidates) == 0) {
  stop(
    "No SAL layer found in ", src_dir, ".\n",
    "Download the Small Area Layer GIS file from DataFirst (login needed):\n  ",
    catalogue_url, "/get-microdata\n",
    "and save it in ", src_dir, " (see the header of this script).",
    call. = FALSE
  )
}

# Prefer a layer whose name mentions SAL / small area
sal_hits <- grepl("sal|small.?area", basename(candidates), ignore.case = TRUE)
if (any(sal_hits)) candidates <- candidates[sal_hits]
if (length(candidates) > 1) {
  stop("More than one possible SAL layer found. Keep one, or set src_layer:\n  ",
       paste(candidates, collapse = "\n  "), call. = FALSE)
}
src_layer <- candidates[1]

if (grepl("\\.gdb$", src_layer, ignore.case = TRUE)) {
  lyr <- st_layers(src_layer)$name
  lyr <- if (length(lyr) > 1) grep("sal", lyr, value = TRUE, ignore.case = TRUE)[1] else lyr
  sal <- st_read(src_layer, layer = lyr, quiet = TRUE)
} else {
  sal <- st_read(src_layer, quiet = TRUE)
}
message("Read ", nrow(sal), " SALs from ", src_layer)


# ============================================================================
# 3. Clip to the ROI and standardise ####
# ============================================================================

if (!exists("cra_roi")) source("R/functions/pipeline_helpers.R")
roi <- cra_roi(catchment_shp, district_shp)

if (is.na(st_crs(sal))) {
  stop("SAL layer has no CRS; check the .prj file in ", src_layer)
}

sal <- sal |>
  st_make_valid() |>
  st_transform(4326)
sal <- sal[lengths(st_intersects(sal, st_as_sfc(roi))) > 0, ]

code_col <- grep("^sal_?code$|^sal$|^sal_?id$", names(sal),
                 value = TRUE, ignore.case = TRUE)[1]
if (is.na(code_col)) {
  stop("No SAL code field found. Fields: ", paste(names(sal), collapse = ", "))
}
sal$sal_code <- as.character(sal[[code_col]])
sal <- sal[, c("sal_code", setdiff(names(sal), "sal_code"))]

if (anyDuplicated(sal$sal_code)) stop("Duplicated SAL codes in ROI subset")

st_write(sal, out, layer = "census2011_sal", delete_dsn = TRUE, quiet = TRUE)


# ============================================================================
# 4. Provenance ####
# ============================================================================

src_files <- if (length(zips)) zips else src_layer
write.csv(
  data.frame(
    source        = "DataFirst: South African Census Community Profiles 2011 (zaf-statssa-saccp-2011-v1)",
    catalogue_url = catalogue_url,
    source_file   = paste(basename(src_files), collapse = "; "),
    md5           = paste(unname(tools::md5sum(src_files[file.exists(src_files) &
                                                   !dir.exists(src_files)])),
                          collapse = "; "),
    layer         = basename(src_layer),
    code_field    = code_col,
    n_sal_roi     = nrow(sal),
    processed_at  = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  file.path(out_dir, "census2011_sal_manifest.csv"),
  row.names = FALSE
)
writeLines(c("[InternetShortcut]", paste0("URL=", catalogue_url)),
           file.path(out_dir, "census2011_sal.url"))

message("Written: ", out, " (", nrow(sal), " SALs)")
