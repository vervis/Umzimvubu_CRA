# ============================================================================
# R/functions/livestock_functions.R
#
# Vulnerability 206 "Livestock density": functions used by
# 109_download_vulnerability_livestock_glw4.R,
# 206_prepare_vulnerability_livestock_density.R and _targets.R (targets
# raw_glw4 and livestock_density).
#
# Source: FAO Gridded Livestock of the World v4 (GLW4), reference year 2020,
# D-DA mapset (dasymetric density, head/km2, ~10 km, EPSG:4326).
#
# Constants
#   GLW4_SPECIES  GLW4 species codes (cattle CTL, sheep SHP, goats GTS)
#   GLW4_URL      download URL pattern (%s = species code)
#   GLW4_CATALOG  FAO catalog record, written to the .url provenance file
#   TLU_WEIGHTS   FAO sub-Saharan Africa TLU weights (cattle 0.7, sheep and
#                 goats 0.1; 1 TLU = 250 kg)
#   TLU_PER_LSU   1 LSU (450 kg, as in the DAFF grazing capacity) = 1.8 TLU
#
# Functions
#   download_glw4()             downloads the GeoTIFFs and writes a manifest
#                               (URL, bytes, md5, access date) and .url file
#   prepare_livestock_density() crops to the ROI, sets negative values to NA,
#                               combines species into TLU/km2 (NA only where
#                               all species are NA) and, if a grazing
#                               capacity raster is given, computes stocking
#                               pressure (LSU/ha / grazing capacity) on the
#                               1 km template. Returns the written paths.
#
# Only function and constant definitions here: tar_source("R/functions")
# loads this folder. Needs the terra package (sf for the roi bbox).
# ============================================================================


GLW4_SPECIES <- c(cattle = "CTL", sheep = "SHP", goats = "GTS")
GLW4_URL     <- "https://storage.googleapis.com/fao-gismgr-glw4-2020-data/DATA/GLW4-2020/MAPSET/D-DA/GLW4-2020.D-DA.%s.tif"
GLW4_CATALOG <- "https://data.apps.fao.org/catalog/iso/9d1e149b-d63f-4213-978b-317a8eb42d02"
TLU_WEIGHTS  <- c(cattle = 0.7, sheep = 0.1, goats = 0.1)
TLU_PER_LSU  <- 450 / 250


# ----------------------------------------------------------------------------
# download_glw4(): GLW4 2020 D-DA GeoTIFFs + provenance
# ----------------------------------------------------------------------------

download_glw4 <- function(out_dir = "data/raw/livestock/glw4_2020",
                          species = GLW4_SPECIES, overwrite = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)

  urls  <- sprintf(GLW4_URL, species)
  files <- file.path(out_dir, basename(urls))

  for (i in seq_along(urls)) {
    if (file.exists(files[i]) && !overwrite) next
    tmp <- paste0(files[i], ".part")
    message("Downloading ", basename(urls[i]))
    utils::download.file(urls[i], tmp, mode = "wb", quiet = TRUE)
    r <- try(terra::rast(tmp), silent = TRUE)
    if (inherits(r, "try-error")) stop(basename(urls[i]), " is not a readable GeoTIFF.")
    file.rename(tmp, files[i])
  }

  writeLines(c("[InternetShortcut]", paste0("URL=", GLW4_CATALOG)),
             file.path(out_dir, "glw4_2020_catalog.url"))

  man <- data.frame(
    species  = names(species),
    code     = unname(species),
    file     = basename(files),
    url      = urls,
    bytes    = file.size(files),
    md5      = unname(tools::md5sum(files)),
    units    = "head/km2",
    mapset   = "D-DA",
    ref_year = 2020,
    accessed = format(Sys.Date())
  )
  utils::write.csv(man, file.path(out_dir, "glw4_2020_manifest.csv"), row.names = FALSE)
  files
}


# ----------------------------------------------------------------------------
# prepare_livestock_density(): heads/km2 -> TLU/km2 (+ stocking pressure)
# ----------------------------------------------------------------------------
#   glw_files : GLW4 D-DA tifs (names in the file must contain CTL/SHP/GTS)
#   roi       : sf bbox (EPSG:4326), the `roi` target
#   template  : 1 km CRA template (only for the stocking-pressure companion)
#   grazing   : grazing capacity raster (LSU/ha, 204); NULL skips the companion
#   buffer    : degrees added around the ROI so bilinear resampling has
#               neighbours at the edge (~2 GLW cells)
# Returns c(tlu, species, stocking) paths (stocking omitted if grazing NULL).

prepare_livestock_density <- function(glw_files, roi, template = NULL,
                                      grazing = NULL,
                                      out_tlu, out_species, out_stocking = NULL,
                                      weights = TLU_WEIGHTS,
                                      buffer = 0.2) {
  code <- vapply(GLW4_SPECIES, function(k) {
    f <- grep(paste0("\\.", k, "\\.tif$"), glw_files, value = TRUE)
    if (length(f) != 1) stop("GLW4 file for ", k, " not found (or not unique).")
    f
  }, "")

  e <- terra::ext(roi["xmin"] - buffer, roi["xmax"] + buffer,
                  roi["ymin"] - buffer, roi["ymax"] + buffer)

  sp <- terra::rast(lapply(code, function(f) terra::crop(terra::rast(f)[[1]], e, snap = "out")))
  names(sp) <- names(GLW4_SPECIES)
  sp <- terra::classify(sp, cbind(-Inf, 0, NA), right = FALSE)

  w   <- weights[names(sp)]
  tlu <- sum(sp * w, na.rm = TRUE)
  tlu <- terra::mask(tlu, sum(!is.na(sp)) > 0, maskvalues = FALSE)
  names(tlu) <- "tlu_per_km2"

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  for (f in c(out_tlu, out_species, out_stocking)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  }
  terra::writeRaster(tlu, out_tlu, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(sp, out_species, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  out <- c(out_tlu, out_species)

  if (!is.null(grazing) && !is.null(out_stocking)) {
    tmpl <- terra::rast(template)
    lsu  <- terra::project(tlu, tmpl, method = "bilinear") / TLU_PER_LSU / 100
    gc   <- terra::project(terra::rast(grazing)[[1]], tmpl, method = "bilinear")
    gc   <- terra::classify(gc, cbind(-Inf, 0, NA), right = TRUE)
    sr   <- terra::mask(lsu / gc, tmpl)
    names(sr) <- "stocking_pressure"
    terra::writeRaster(sr, out_stocking, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
    out <- c(out, out_stocking)
  }
  out
}
