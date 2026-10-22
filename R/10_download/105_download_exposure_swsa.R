# ============================================================================
# 105_download_exposure_swsa.R
#
# Exposure 212: Strategic Water Source Areas (SANBI BGIS) for the CRA
# rectangle (cra_roi()).
#   Surface water (hybrid):
#     National: 2021 surface-water SWSA update (service
#       2021SWSAsForSurfaceWater, layer 9 SWSA_sw_March2021_DD, incl. Lesotho
#       and eSwatini), a fine-scale delineation of the 22 National SWSAs.
#     Sub-national: 2017 SWSA surface polygons with Criteria Sub-national
#       (here *Pondoland Coast). These high-runoff coastal catchments were left
#       out of the National set as "water secure", and the 2021 update only
#       re-delineated the National set.
#     prepare_swsa() takes the maximum weight where they overlap.
#   Groundwater: 2017 SWSA (Le Maitre et al. 2018; service 2017SWSA, layer 0),
#     Type gw + swgw. There is no 2021 groundwater update.
# Both are written to data/raw/swsa/Strategic_Water_Source_Areas_CRA.gpkg
# (layers SWSA_surface_water, SWSA_groundwater) with the fields prepare_swsa()
# uses: Criteria (National / Sub-national), Name, Type, source.
# Requests send a browser user agent: the BGIS server returns 403 to the
# httr2 default.
# ============================================================================

library(sf)
library(httr2)
library(dplyr)

source("R/functions/pipeline_helpers.R")

roi  <- cra_roi()
gpkg <- "data/raw/swsa/Strategic_Water_Source_Areas_CRA.gpkg"
dir.create(dirname(gpkg), recursive = TRUE, showWarnings = FALSE)

bgis <- "https://bgismaps.sanbi.org/server/rest/services/BGIS_Projects/"
url_sw <- paste0(bgis, "2021SWSAsForSurfaceWater/MapServer/9/query")
url_gw <- paste0(bgis, "2017SWSA/MapServer/0/query")

get_swsa <- function(url) {
  r <- request(url) |>
    req_user_agent("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/141.0.0.0 Safari/537.36") |>
    req_headers(Accept = "application/json,text/html;q=0.9,*/*;q=0.8",
                `Accept-Language` = "en-ZA,en;q=0.9") |>
    req_url_query(
      where = "1=1",
      geometry = paste(roi[c("xmin", "ymin", "xmax", "ymax")], collapse = ","),
      geometryType = "esriGeometryEnvelope", inSR = 4326,
      spatialRel = "esriSpatialRelIntersects",
      outFields = "*", returnGeometry = "true", outSR = 4326, f = "geojson"
    ) |>
    req_retry(max_tries = 5) |>
    req_perform()
  tmp <- tempfile(fileext = ".geojson")
  writeBin(resp_body_raw(r), tmp)
  st_crop(st_make_valid(st_read(tmp, quiet = TRUE)), roi)
}

s17 <- get_swsa(url_gw)

sw <- rbind(
  get_swsa(url_sw) |>
    transmute(Criteria = Priority, Name = SWSAName, Type = "sw",
              source = "SANBI BGIS 2021 surface-water SWSA (March 2021)"),
  s17 |>
    filter(Type %in% c("sw", "swgw"), Criteria == "Sub-national") |>
    transmute(Criteria, Name = Name_sw, Type = "sw",
              source = "SANBI BGIS 2017 SWSA, sub-national surface water")
)
gw <- s17 |>
  filter(Type %in% c("gw", "swgw")) |>
  transmute(Criteria, Name = Name_gw, Type,
            source = "SANBI BGIS 2017 SWSA (Le Maitre et al. 2018)")

if (file.exists(gpkg)) file.remove(gpkg)
st_write(sw, gpkg, layer = "SWSA_surface_water", quiet = TRUE)
st_write(gw, gpkg, layer = "SWSA_groundwater", append = TRUE, quiet = TRUE)

write.csv(data.frame(layer = c("SWSA_surface_water", "SWSA_groundwater"),
                     source = c(paste(url_sw, "+ Sub-national from", url_gw), url_gw),
                     n = c(nrow(sw), nrow(gw)),
                     accessed = format(Sys.Date())),
          file.path(dirname(gpkg), "swsa_manifest.csv"), row.names = FALSE)

print(st_drop_geometry(sw))
print(st_drop_geometry(gw))
