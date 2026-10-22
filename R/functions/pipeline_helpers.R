# ============================================================================
# R/functions/pipeline_helpers.R
#
# Shared constants and helpers for _targets.R and the stage scripts. Only
# function and constant definitions: tar_source("R/functions") loads them.
#
# CRA extent: one rectangle for every layer, the bounding box of the
#   Umzimvubu catchment + Alfred Nzo district in SA Albers Equal Area,
#   snapped outwards to whole 1 km cells. cra_extent() returns it in the CRA
#   CRS; cra_roi() returns the WGS84 bbox that encloses it (for downloads and
#   crops of geographic sources), so every source covers the full template.
#
# Constants
#   CRA_CRS         SA Albers Equal Area (SANLC 2022 CRS)
#   CRA_RES         template cell size (m)
#   CRA_BOUNDARIES  catchment and district shapefiles
#
# Functions
#   cra_extent()  terra SpatExtent of the CRA rectangle (CRA_CRS)
#   cra_roi()     sf bbox (EPSG:4326) enclosing the CRA rectangle
#   run_script()  runs a stage script inside a target and checks its outputs
# ============================================================================

CRA_CRS <- paste(
  "+proj=aea +lat_0=0 +lon_0=25 +lat_1=-24 +lat_2=-33",
  "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)
CRA_RES <- 1000
CRA_BOUNDARIES <- c(
  catchment = "data/raw/boundaries/Umzimvubu_Catchment.shp",
  district  = "data/raw/boundaries/Alfred_Nzo.shp"
)


cra_extent <- function(catchment_shp = CRA_BOUNDARIES[["catchment"]],
                       district_shp = CRA_BOUNDARIES[["district"]],
                       crs = CRA_CRS, res = CRA_RES) {
  e <- sapply(c(catchment_shp, district_shp), function(f) {
    as.vector(terra::ext(terra::project(terra::vect(f), crs)))
  })
  terra::ext(floor(min(e[1, ]) / res) * res, ceiling(max(e[2, ]) / res) * res,
             floor(min(e[3, ]) / res) * res, ceiling(max(e[4, ]) / res) * res)
}


cra_roi <- function(catchment_shp = CRA_BOUNDARIES[["catchment"]],
                    district_shp = CRA_BOUNDARIES[["district"]],
                    crs = CRA_CRS, res = CRA_RES) {
  e <- as.vector(cra_extent(catchment_shp, district_shp, crs, res))
  p <- sf::st_as_sfc(sf::st_bbox(c(xmin = e[[1]], xmax = e[[2]],
                                   ymin = e[[3]], ymax = e[[4]]),
                                 crs = sf::st_crs(crs)))
  p <- sf::st_segmentize(p, res)
  sf::st_bbox(sf::st_transform(p, 4326))
}


run_script <- function(script, outputs, ...) {
  force(list(...))
  old_opts <- options(browser = function(...) invisible(NULL))
  on.exit(options(old_opts), add = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  message("Running ", script)
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env, chdir = FALSE)

  missing <- outputs[!file.exists(outputs)]
  if (length(missing)) {
    stop(basename(script), " finished but did not create:\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
  outputs
}
