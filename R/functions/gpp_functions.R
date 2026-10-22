# ============================================================================
# R/functions/gpp_functions.R
#
# Vulnerability 208 "Biomass availability" (sensitivity alternative to
# grazing_capacity): Global Pasture Watch pantropical calibrated grassland
# GPP, 30 m, annual, gC m-2 yr-1 (UInt16), v1.2 (Xu, Isik, Parente et al.;
# Zenodo, CC BY 4.0). uGPP x land-cover-specific light use efficiency from
# the Global Harmonized Land Cover map of the same year.
#
# Constants
#   GPP_RECORDS  Zenodo record id and file name per year (2022-2024)
#
# Functions
#   download_gpp()        windowed read of each global COG over the CRA
#                         rectangle (GDAL /vsicurl/, HTTP range requests, so
#                         only the ROI tiles are transferred); writes one
#                         GeoTIFF per year and a manifest
#   prepare_biomass_gpp() mean annual grassland GPP across years on the 30 m
#                         grid (0 = no grassland -> NA); stage 30 averages it
#                         to the template
# ============================================================================


GPP_RECORDS <- data.frame(
  year   = 2022:2024,
  record = c(23063228, 23063230, 23063232),
  file   = sprintf("gpw_gpp.grassland_lue.model_m_30m_s_%d0101_%d1231_go03_epsg.4326_v1.2.tif",
                   2022:2024, 2022:2024)
)


download_gpp <- function(roi, out_dir = "data/raw/gpp", records = GPP_RECORDS,
                         overwrite = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # the Zenodo URLs end in /content, so no extension filter
  ext_old <- Sys.getenv("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", NA)
  Sys.unsetenv("CPL_VSIL_CURL_ALLOWED_EXTENSIONS")
  if (!is.na(ext_old)) on.exit(Sys.setenv(CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ext_old), add = TRUE)
  Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR", GDAL_HTTP_MAX_RETRY = "5",
             GDAL_HTTP_RETRY_DELAY = "5")

  e <- terra::ext(roi[["xmin"]], roi[["xmax"]], roi[["ymin"]], roi[["ymax"]])
  records$url <- sprintf("https://zenodo.org/api/records/%d/files/%s/content",
                         records$record, records$file)
  records$out <- file.path(out_dir, sprintf("gpw_grassland_gpp_%d_cra.tif", records$year))

  for (i in seq_len(nrow(records))) {
    if (file.exists(records$out[i]) && !overwrite) next
    message("GPP ", records$year[i], ": reading ROI window from Zenodo")
    r <- terra::rast(paste0("/vsicurl/", records$url[i]))
    r <- terra::crop(r, e, snap = "out")
    names(r) <- paste0("gpp_", records$year[i])
    terra::writeRaster(r, records$out[i], overwrite = TRUE, datatype = "INT2U",
                       gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  }

  utils::write.csv(data.frame(
    year = records$year, zenodo_record = records$record, url = records$url,
    file = records$out, bytes = file.size(records$out),
    md5 = unname(tools::md5sum(records$out)), units = "gC m-2 yr-1",
    licence = "CC BY 4.0", accessed = format(Sys.Date())
  ), file.path(out_dir, "gpp_manifest.csv"), row.names = FALSE)
  records$out
}


prepare_biomass_gpp <- function(gpp_files, out) {
  r <- terra::rast(gpp_files)
  r <- terra::classify(r, cbind(0, NA))
  m <- terra::mean(r, na.rm = TRUE)
  names(m) <- "grassland_gpp_mean"
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(m, out, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  out
}
