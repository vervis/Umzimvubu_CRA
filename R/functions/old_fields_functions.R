# ============================================================================
# R/functions/old_fields_functions.R
#
# Vulnerability 207 "Old fields": land in lower ecological condition. Fallow
# land and old fields (SANLC 2022 classes 42-46 and 73) are former cropland
# with degraded soils and vegetation cover, so they are more sensitive to
# drought (lower infiltration and water retention, slower recovery, more
# erosion). Indicator = share of the 1 km cell under old fields.
#
# prepare_old_fields(): 20 m 0/1 old-field mask averaged to the template;
# also writes the 20 m old-field classes (0 elsewhere) for mapping. Returns
# c(share, classes).
# ============================================================================


SANLC_OLD_FIELD_CLASSES <- c(42:46, 73)


prepare_old_fields <- function(sanlc_file, template, out_share, out_classes,
                               of_classes = SANLC_OLD_FIELD_CLASSES) {
  s <- terra::rast(sanlc_file)[[1]]
  s <- terra::classify(s, cbind(0, NA))

  cls <- terra::ifel(s %in% of_classes, s, 0)
  names(cls) <- "old_field_class"

  tmpl  <- terra::rast(template)
  share <- terra::project(terra::ifel(s %in% of_classes, 1, 0), tmpl, method = "average")
  share <- terra::mask(share, tmpl)
  names(share) <- "old_fields_share"

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  for (f in c(out_share, out_classes)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(share, out_share, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(cls, out_classes, overwrite = TRUE, datatype = "INT1U", gdal = gdal)
  c(out_share, out_classes)
}
