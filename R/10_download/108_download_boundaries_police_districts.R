# ============================================================================
# 108_download_boundaries_police_districts.R
#
# Stats SA police-district boundaries, the geography of the Census 2022
# police-district tables (DataFirst casd-2022-policedist v1). Used by 205
# (unimproved water, Census 2022).
#
# MANUAL STEP: downloaded by hand from the Stats SA website (30 Sep 2026)
# and saved, unchanged, as
#   data/raw/boundaries/police_districts/PoliceDistrict.zip
# Contents: PoliceDistrict.shp (EPSG:4326), 1,167 districts; SAPS component
# boundaries CREATE_DT 20241001, VERSION 1.3.5; fields PD_CODE, PD_Name
# (= COMPNT_NM), MN_MDB_C / MN_NAME (municipality), DC_MDB_C.
# PD_Name matches the police_district names of the DataFirst tables one to
# one after police_key(), so no name lookup is needed.
#
# This script does not download anything. It checks the zip, writes a
# provenance manifest and lists the districts that touch the CRA rectangle.
# prepare_unimproved_water_2022() reads the zip directly.
#
# Outputs (data/raw/boundaries/police_districts/):
#   police_districts_manifest.csv   file, checksum, feature count, ROI districts
#   police_districts.url            Stats SA link (edit if you have the exact page)
#
# _targets.R entry (20 Prepare section, like the Census 2011 zip):
#   tar_target(raw_police_districts,
#              "data/raw/boundaries/police_districts/PoliceDistrict.zip",
#              format = "file"),
# ============================================================================

library(sf)

police_zip    <- "data/raw/boundaries/police_districts/PoliceDistrict.zip"
out_dir       <- dirname(police_zip)
catchment_shp <- "data/raw/boundaries/Umzimvubu_Catchment.shp"
district_shp  <- "data/raw/boundaries/Alfred_Nzo.shp"
source_url    <- "https://www.statssa.gov.za/?p=18624"  

stopifnot(file.exists(police_zip))
if (!exists("cra_roi")) source("R/functions/pipeline_helpers.R")
if (!exists("read_police_districts")) source("R/functions/water_access_functions.R")

roi <- cra_roi(catchment_shp, district_shp)
all_pd <- read_police_districts(police_zip)
roi_pd <- read_police_districts(police_zip, roi)
roi_pd <- roi_pd[lengths(st_intersects(roi_pd, st_as_sfc(roi))) > 0, ]

write.csv(
  data.frame(
    source        = "Stats SA police-district boundaries (SAPS components, CREATE_DT 20241001, v1.3.5), downloaded by hand",
    source_url    = source_url,
    file          = basename(police_zip),
    md5           = unname(tools::md5sum(police_zip)),
    n_districts   = nrow(all_pd),
    n_in_roi      = nrow(roi_pd),
    roi_districts = paste(sort(roi_pd$pd_name), collapse = "; "),
    checked_at    = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  file.path(out_dir, "police_districts_manifest.csv"),
  row.names = FALSE
)
writeLines(c("[InternetShortcut]", paste0("URL=", source_url)),
           file.path(out_dir, "police_districts.url"))

message(nrow(all_pd), " police districts; ", nrow(roi_pd), " touch the ROI:")
print(sort(roi_pd$pd_name))
