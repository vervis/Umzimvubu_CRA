# ============================================================================
# R/functions/swsa_functions.R
#
# Exposure 212 "Strategic Water Source Areas" and the SWSA degradation
# layers. Used by 212_prepare_exposure_swsa.R and _targets.R (targets swsa
# and swsa_degradation).
#
# Sources (SANBI BGIS, downloaded by 105 into
# data/raw/swsa/Strategic_Water_Source_Areas_CRA.gpkg):
#   SWSA_surface_water  2021 surface-water SWSA update (National) + 2017
#                       Sub-national surface SWSAs (Pondoland Coast)
#   SWSA_groundwater    2017 SWSA groundwater, Type gw + swgw (Le Maitre et
#                       al. 2018, WRC TT 743/1/18; National / Sub-national)
#
# Framing: SWSA = exposed natural asset (regional water supply), higher =
# more exposure. Condition of the asset (alien trees, old fields, grazing)
# is vulnerability, so it is NOT built into 212. The degradation layers put
# the two together for NbS targeting, next to the composite, not inside it.
#
# Constants
#   SWSA_LAYERS        gpkg layer names (surface, ground)
#   SWSA_TIER_WEIGHTS  weight per Criteria value (National 1, Sub-national 0.5)
#
# Functions
#   swsa_layer_share()        one layer -> tier-weighted share of each 1 km cell
#                             (rasterised at res/fine, averaged to template)
#   prepare_swsa()            surface + ground shares and the combined
#                             indicator E = w_sw * s_sw + (1 - w_sw) * s_gw.
#                             Returns c(combined, surface, ground).
#   prepare_swsa_degradation() s_sw x alien tree cover and s_gw x land-use
#                             pressure (default old fields share).
#                             Returns c(surface, ground).
# ============================================================================


SWSA_LAYERS <- c(surface = "SWSA_surface_water", ground = "SWSA_groundwater")
SWSA_TIER_WEIGHTS <- c("National" = 1, "Sub-national" = 0.5)


swsa_layer_share <- function(gpkg, layer, template, tier_w = SWSA_TIER_WEIGHTS,
                             fine = 10) {
  tmpl <- terra::rast(template)
  v <- terra::makeValid(terra::project(terra::vect(gpkg, layer = layer), tmpl))
  miss <- setdiff(unique(v$Criteria), names(tier_w))
  if (length(miss)) stop(layer, ": no tier weight for Criteria ",
                         paste(miss, collapse = ", "))
  v$w <- unname(tier_w[v$Criteria])
  f <- terra::disagg(terra::rast(tmpl), fine)
  r <- terra::rasterize(v, f, field = "w", fun = "max", background = 0)
  r <- terra::aggregate(r, fine, fun = "mean")
  terra::mask(r, tmpl)
}


prepare_swsa <- function(gpkg, template, out_combined, out_surface, out_ground,
                         w_sw = 0.5, tier_w = SWSA_TIER_WEIGHTS, fine = 10) {
  sw <- swsa_layer_share(gpkg, SWSA_LAYERS[["surface"]], template, tier_w, fine)
  gw <- swsa_layer_share(gpkg, SWSA_LAYERS[["ground"]], template, tier_w, fine)
  e  <- w_sw * sw + (1 - w_sw) * gw
  names(sw) <- "swsa_surface"
  names(gw) <- "swsa_ground"
  names(e)  <- "swsa_combined"

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  out <- c(combined = out_combined, surface = out_surface, ground = out_ground)
  for (f in out) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(e,  out[["combined"]], overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(sw, out[["surface"]],  overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(gw, out[["ground"]],   overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  out
}


# Product of SWSA importance (0-1) and pressure (share of cell, 0-1).
# Cells outside a SWSA are 0 even where the pressure layer is NA; inside a
# SWSA, NA pressure (e.g. Alfred Nzo outside the MapWAPS extent) stays NA.
prepare_swsa_degradation <- function(surface, ground, alien, gw_pressure,
                                     template, out_surface, out_ground) {
  tmpl <- terra::rast(template)
  on_tmpl <- function(f) {
    r <- terra::rast(f)[[1]]
    if (!terra::compareGeom(r, tmpl, stopOnError = FALSE)) {
      r <- terra::project(r, tmpl, method = "average")
    }
    r
  }
  deg <- function(s, p) {
    s <- on_tmpl(s)
    terra::mask(terra::ifel(s == 0, 0, s * on_tmpl(p)), tmpl)
  }
  d_sw <- deg(surface, alien)
  d_gw <- deg(ground, gw_pressure)
  names(d_sw) <- "swsa_surface_alien_pressure"
  names(d_gw) <- "swsa_ground_landuse_pressure"

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  out <- c(surface = out_surface, ground = out_ground)
  for (f in out) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(d_sw, out[["surface"]], overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  terra::writeRaster(d_gw, out[["ground"]],  overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  out
}
