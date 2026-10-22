# ============================================================================
# 106_download_boundaries_mdb_wards.R
#
# Ward boundaries from the Municipal Demarcation Board (MDB), downloaded in
# one step from the MDB's public ArcGIS Online feature services and clipped
# to the CRA region of interest (Umzimvubu catchment + Alfred Nzo bounding
# box).
#
# Why this source: the MDB ward-delimitation web page offers one file per
# municipality or province. The MDB also publishes every ward vintage as a
# national feature service, so one query returns all wards in the ROI.
#
#   2020 delimitation: used for the 2021 local government elections AND for
#     Census 2022 ward-level outputs. Needed for option 2 of the rainfed
#     crop reliance indicator (smallholder cropland per household).
#     The MDB item is titled "MDB Wards 2021"; the service is MDB_Wards_2020.
#   2011 delimitation: matches Census 2011 wards (SALs nest within them).
#     Optional. Useful for summarising option 3 (SAL level) to wards.
#     Known issue: the MDB_Wards_2011 service holds every ward twice (8,554
#     features for 4,277 wards; e.g. WardID 21001001 is OBJECTID 773 and
#     5050, with identical attributes and area). The script checks that the
#     copies are identical and keeps one. Its municipality names/codes are
#     post-2016 ones (e.g. EC101 = Dr Beyers Naude, not Camdeboo), so join
#     Census 2011 data on WardID, not on municipality name.
#   2016 delimitation: used for the 2016 local government elections and for
#     Community Survey 2016 ward tables on "Geography 2016" in SuperWEB2
#     (dataset Householdagriculture16_wd). 4,392 wards nationally, no
#     duplicates. Same field names as the 2011 service.
#   WardIDs look alike across vintages (e.g. 24401001) but the boundaries
#   differ, so always join a table to the ward vintage it was tabulated on.
#
# Outputs (data/raw/boundaries/mdb_wards_<vintage>/):
#   mdb_wards_<vintage>_cra.gpkg     wards intersecting the ROI, EPSG:4326
#   mdb_wards_<vintage>_manifest.csv source, query date, feature count
#   mdb_wards_<vintage>.url          link to the source service
#
# Standardised fields: ward_id, ward_no, lm_code, lm_name, dm_code, dm_name,
# province, plus the original MDB fields.
#
# _targets.R entry (10 Download section):
#   download_target(
#     "raw_wards_2020",
#     "R/10_download/106_download_boundaries_mdb_wards.R",
#     "data/raw/boundaries/mdb_wards_2020/mdb_wards_2020_cra.gpkg"
#   ),
# ============================================================================

library(sf)
library(httr2)

# ============================================================================
# 1. Settings ####
# ============================================================================

# Vintages to download ("2011", "2016", "2020")
vintages <- c("2011", "2016", "2020")

mdb_base <- "https://services7.arcgis.com/oeoyTUJC8HEeYsRB/arcgis/rest/services"

ward_services <- list(
  "2020" = list(
    service = "MDB_Wards_2020",
    oid     = "FID",
    fields  = c(ward_id = "WardID", ward_no = "WardNo",
                lm_code = "CAT_B", lm_name = "Municipali",
                dm_code = "DistrictCo", dm_name = "District",
                province = "Province")
  ),
  "2011" = list(
    service = "MDB_Wards_2011",
    oid     = "OBJECTID",
    fields  = c(ward_id = "WardID", ward_no = "WardNumber",
                lm_code = "LocalMunicipalityCode",
                lm_name = "LocalMunicipalityName",
                dm_code = "DistrictMunicipalityCode",
                dm_name = "DistrictMunicipalityName",
                province = "ProvinceName")
  ),
  "2016" = list(
    service = "MDB_Wards_2016",
    oid     = "OBJECTID",
    fields  = c(ward_id = "WardID", ward_no = "WardNumber",
                lm_code = "LocalMunicipalityCode",
                lm_name = "LocalMunicipalityName",
                dm_code = "DistrictMunicipalityCode",
                dm_name = "DistrictMunicipalityName",
                province = "ProvinceName")
  )
)

catchment_shp <- "data/raw/boundaries/Umzimvubu_Catchment.shp"
district_shp  <- "data/raw/boundaries/Alfred_Nzo.shp"

page_size <- 1000   # service maxRecordCount is 2000


# ============================================================================
# 2. ROI (same bounding box as the rest of the pipeline) ####
# ============================================================================

if (!exists("cra_roi")) source("R/functions/pipeline_helpers.R")
roi <- cra_roi(catchment_shp, district_shp)

envelope <- paste(roi[c("xmin", "ymin", "xmax", "ymax")], collapse = ",")


# ============================================================================
# 3. Query helpers ####
# ============================================================================

query_layer <- function(service, ...) {
  request(sprintf("%s/%s/FeatureServer/0/query", mdb_base, service)) |>
    req_url_query(
      where        = "1=1",
      geometry     = envelope,
      geometryType = "esriGeometryEnvelope",
      inSR         = 4326,
      spatialRel   = "esriSpatialRelIntersects",
      ...
    ) |>
    req_retry(max_tries = 5, backoff = function(i) 2^i) |>
    req_timeout(120)
}

count_features <- function(service) {
  resp <- query_layer(service, returnCountOnly = "true", f = "json") |>
    req_perform() |>
    resp_body_json()
  if (!is.null(resp$error)) stop(service, ": ", resp$error$message)
  resp$count
}

# Pages through the service; returns one sf object in EPSG:4326.
# Pages are ordered by the object ID, which is unique, so no feature is
# skipped or repeated between pages (WardID is not unique in the 2011 service).
download_features <- function(service, n_expected, oid) {
  pages  <- list()
  offset <- 0
  repeat {
    txt <- query_layer(
      service,
      outFields         = "*",
      returnGeometry    = "true",
      outSR             = 4326,
      orderByFields     = oid,
      resultOffset      = offset,
      resultRecordCount = page_size,
      f                 = "geojson"
    ) |>
      req_perform() |>
      resp_body_string()

    page <- st_read(txt, quiet = TRUE)
    if (nrow(page) == 0) break
    pages[[length(pages) + 1]] <- page
    offset <- offset + nrow(page)
    if (offset >= n_expected) break
  }
  do.call(rbind, pages)
}

# Adds standard names alongside the original fields (partial matching guards
# against truncated shapefile-style names such as "Municipali")
standardise_fields <- function(x, map) {
  for (new in names(map)) {
    hit <- grep(paste0("^", map[[new]]), names(x), value = TRUE,
                ignore.case = TRUE)[1]
    if (is.na(hit)) {
      warning("Field '", map[[new]], "' not found; ", new, " set to NA")
      x[[new]] <- NA_character_
    } else {
      x[[new]] <- as.character(x[[hit]])
    }
  }
  x[, c(names(map), setdiff(names(x), names(map)))]
}

# GeoPackage field names are case-insensitive, and "fid" is reserved for the
# GeoPackage's own feature ID. MDB_Wards_2020 has "Province" (clashes with
# the standard "province") and "FID" (clashes with fid), which makes
# st_write() fail with "Field count reached: duplicate names present?".
# Original fields that clash are renamed with an "mdb_" prefix.
gpkg_safe_names <- function(x, standard) {
  geom <- attr(x, "sf_column")
  nms  <- names(x)
  orig <- !nms %in% c(standard, geom)
  clash <- orig & tolower(nms) %in% c(tolower(standard), "fid")
  names(x)[clash] <- paste0("mdb_", nms[clash])
  if (anyDuplicated(tolower(names(x)))) {
    stop("Field names still clash (case-insensitive): ",
         paste(names(x)[duplicated(tolower(names(x)))], collapse = ", "))
  }
  x
}

# Drops repeated copies of the same ward. Stops if the copies differ in any
# attribute (other than the object ID) or in area, because then they are not
# true duplicates and need a closer look.
drop_duplicate_wards <- function(x, service, oid) {
  dup_ids <- unique(x$ward_id[duplicated(x$ward_id)])
  if (length(dup_ids) == 0) return(x)

  attrs <- setdiff(names(st_drop_geometry(x)), oid)
  dups  <- x[x$ward_id %in% dup_ids, ]
  dups$area_m2 <- round(as.numeric(st_area(dups)))

  n_versions <- tapply(
    seq_len(nrow(dups)), dups$ward_id,
    function(i) nrow(unique(st_drop_geometry(dups[i, c(attrs, "area_m2")])))
  )
  if (any(n_versions > 1)) {
    stop(service, ": ", sum(n_versions > 1), " ward_id values have copies ",
         "that differ, e.g. ", names(n_versions)[n_versions > 1][1])
  }

  message(service, ": ", length(dup_ids), " wards appear more than once ",
          "(identical copies); keeping one copy of each")
  x[!duplicated(x$ward_id), ]
}


# ============================================================================
# 4. Download each vintage ####
# ============================================================================

for (v in vintages) {
  svc     <- ward_services[[v]]
  out_dir <- file.path("data/raw/boundaries", paste0("mdb_wards_", v))
  out     <- file.path(out_dir, sprintf("mdb_wards_%s_cra.gpkg", v))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n <- count_features(svc$service)
  message(svc$service, ": ", n, " wards intersect the ROI")

  wards <- download_features(svc$service, n, svc$oid) |>
    st_make_valid() |>
    standardise_fields(svc$fields)

  # n counts features, not wards, so check completeness before de-duplicating
  if (nrow(wards) != n) {
    stop(svc$service, ": expected ", n, " features, downloaded ", nrow(wards))
  }
  wards <- drop_duplicate_wards(wards, svc$service, svc$oid)
  message(svc$service, ": ", nrow(wards), " unique wards")

  if (anyDuplicated(wards$ward_id)) {
    stop(svc$service, ": duplicated ward_id values")
  }

  wards <- gpkg_safe_names(wards, names(svc$fields))
  st_write(wards, out, layer = paste0("mdb_wards_", v),
           delete_dsn = TRUE, delete_layer = TRUE, quiet = TRUE)

  # Provenance (tracked by git: *_manifest.csv and *.url)
  write.csv(
    data.frame(
      source        = "Municipal Demarcation Board ArcGIS Online",
      service_url   = sprintf("%s/%s/FeatureServer/0", mdb_base, svc$service),
      ward_vintage  = v,
      roi_bbox_4326 = envelope,
      n_features    = n,
      n_wards       = nrow(wards),
      n_lm          = length(unique(wards$lm_code)),
      downloaded_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    ),
    file.path(out_dir, sprintf("mdb_wards_%s_manifest.csv", v)),
    row.names = FALSE
  )
  writeLines(
    c("[InternetShortcut]",
      sprintf("URL=%s/%s/FeatureServer/0", mdb_base, svc$service)),
    file.path(out_dir, sprintf("mdb_wards_%s.url", v))
  )

  # QA
  print(table(wards$lm_code, useNA = "ifany"))
  message("Written: ", out)
}
