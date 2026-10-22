# ============================================================================
# 000_create_directories.R
#
# Creates the project folder structure: stage script folders, raw data
# folders per source and the processed / output folders the pipeline writes
# to. Existing folders are left unchanged; prints which ones were created.
# ============================================================================

dirs <- c(
  file.path("R", c("00_setup", "10_download", "20_prepare", "30_indicators",
                   "40_components", "50_risk", "60_validation",
                   "70_visualisation", "80_export", "90_outputs", "functions")),
  file.path("data/raw", c("boundaries", "nex_gddp_cmip6", "landcover/sanlc2022",
                          "grazing_capacity", "meta_rwi", "swsa", "census2011",
                          "census2022", "dws", "livestock", "mapwaps", "rivers", "gpp")),
  file.path("data/processed", c("template", "landcover/sanlc2022", "hazard",
                                "vulnerability", "exposure", "nbs", "indicators",
                                "components", "risk")),
  "outputs/tables", "outputs/maps"
)

created <- vapply(dirs, dir.create, logical(1), recursive = TRUE, showWarnings = FALSE)
print(data.frame(directory = dirs, created_now = created, exists = dir.exists(dirs)))
