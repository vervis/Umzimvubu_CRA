# ============================================================================
# _targets.R: Umzimvubu Climate Risk Assessment pipeline
#
# Stages follow the R/ folder numbering:
#   00 inputs > 10 download > 20 prepare > 30 indicators > 40 components >
#   50 risk > 60 validation > 70 visualisation > 80 export
#
# Every layer uses one extent: the CRA rectangle (cra_extent() / cra_roi() in
# R/functions/pipeline_helpers.R), the bounding box of the Umzimvubu
# catchment + Alfred Nzo district snapped to the 1 km template.
#
# Everyday use (from the project root):
#   targets::tar_visnetwork()   # see the pipeline and what is out of date
#   targets::tar_make()         # run only the steps that are out of date
#   targets::tar_read(roi)      # load a target's value
#
# Downloads are not re-run unless RUN_DOWNLOADS is TRUE: with FALSE the raw
# files are tracked, so run the 10_download scripts by hand first.
# ============================================================================

library(targets)
library(tarchetypes)

tar_option_set(
  packages = c("sf", "terra", "dplyr", "purrr", "stringr", "readr", "haven", "readxl"),
  format   = "rds"
)

tar_source("R/functions")

RUN_DOWNLOADS <- FALSE

download_target <- function(name, script, outputs) {
  command <- if (RUN_DOWNLOADS) {
    bquote(run_script(.(script), .(outputs)))
  } else {
    as.call(c(quote(c), as.list(outputs)))
  }
  tar_target_raw(name, command, format = "file")
}

script_target <- function(name, path) {
  tar_target_raw(name, as.call(list(quote(c), path)), format = "file")
}

# NEX-GDDP periods built in stage 30; each must exist as "{period}" in the
# file names written by 201. Keep NEX_SSPS equal to use_ssps in 201.
NEX_SSPS <- c("ssp126", "ssp245")
NEX_PERIODS <- c(
  "hist_1981_2014",
  as.vector(outer(c("near_2021_2040_", "mid_2041_2070_", "late_2071_2100_"),
                  NEX_SSPS, paste0))
)

INDICATOR_REGISTER <- "R/30_indicators/indicator_register.csv"
indicator_register <- read_indicator_register(INDICATOR_REGISTER)

indicator_targets <- tar_map(
  values = indicator_register,
  names  = "id",
  tar_target(
    indicator,
    make_indicator(
      id = id, component = component, src = src,
      temporal = temporal, periods = NEX_PERIODS,
      template = template_grid, resample = resample,
      direction = direction, normalise = normalise, trim = trim,
      lower = lower, upper = upper, na_fill = na_fill,
      transform = transform,
      upstream = upstream
    ),
    format = "file"
  )
)

SANLC <- "data/processed/landcover/sanlc2022/sanlc2022_umzimvubu_alfred_nzo_extent.tif"
VUL   <- "data/processed/vulnerability"
EXP   <- "data/processed/exposure"


list(

  # ==========================================================================
  # 00 Inputs
  # ==========================================================================
  tar_target(catchment_shp, CRA_BOUNDARIES[["catchment"]], format = "file"),
  tar_target(district_shp, CRA_BOUNDARIES[["district"]], format = "file"),
  tar_target(roi, cra_roi(catchment_shp, district_shp)),
  tar_target(
    template_grid,
    make_template_grid(catchment_shp, district_shp,
                       out = "data/processed/template/cra_template_1km.tif"),
    format = "file"
  ),

  # ==========================================================================
  # 10 Download
  # ==========================================================================
  download_target("raw_nex_gddp",
                  "R/10_download/101_download_hazards_all_ssps_common_ensemble.R",
                  "data/raw/nex_gddp_cmip6"),
  download_target("sanlc_cra",
                  "R/10_download/102_download_vulnerability_agriculture_layers.R",
                  SANLC),
  download_target("raw_grazing",
                  "R/10_download/103_download_vulnerability_grazing_capacity.R",
                  "data/raw/grazing_capacity/grazing_capacity_2016_sa.gpkg"),
  download_target("rwi",
                  "R/10_download/104_download_vulnerability_meta_rwi.R",
                  file.path(VUL, "meta_rwi/meta_relative_wealth_index_roi.tif")),
  download_target("raw_swsa",
                  "R/10_download/105_download_exposure_swsa.R",
                  "data/raw/swsa/Strategic_Water_Source_Areas_CRA.gpkg"),
  download_target("raw_glw4",
                  "R/10_download/109_download_vulnerability_livestock_glw4.R",
                  file.path("data/raw/livestock/glw4_2020",
                            sprintf("GLW4-2020.D-DA.%s.tif", c("CTL", "SHP", "GTS")))),
  download_target("raw_mapwaps",
                  "R/10_download/110_download_vulnerability_mapwaps.R",
                  c("data/raw/mapwaps/MAPWAPS_uMzimvubu_19082025.tif",
                    "data/raw/mapwaps/mapwaps_class_lookup.csv")),
  download_target("raw_rivers",
                  "R/10_download/111_download_vulnerability_rivers.R",
                  "data/raw/rivers/nfepa_rivers_cra.gpkg"),
  download_target("raw_gpp",
                  "R/10_download/112_download_vulnerability_gpp.R",
                  sprintf("data/raw/gpp/gpw_grassland_gpp_%d_cra.tif", 2022:2024)),

  # Census and boundary files downloaded by hand
  tar_target(raw_census2011_sal_zip, "data/raw/census2011/mapdata-sa.zip", format = "file"),
  tar_target(raw_census2011_water_source, "data/raw/census2011/source-of-water.dta",
             format = "file"),
  tar_target(raw_census2011_piped_water, "data/raw/census2011/piped-water.dta",
             format = "file"),
  tar_target(raw_census2011_income, "data/raw/census2011/annual-household-income.dta",
             format = "file"),
  tar_target(raw_census2011_activity_union,
             "data/raw/census2011/Census2011_Ward2016_ActivityUnion_HH.csv", format = "file"),
  tar_target(raw_census2011_agric_sample,
             "data/raw/census2011/sa-census-2011-agricultural-households-v1-20151104-stata.zip",
             format = "file"),
  tar_target(raw_census2022_water_source,
             "data/raw/census2022/police-district-v1/Police-District/casd-2022-policedist-stata/c2022-watersource-policedist-v1.dta",
             format = "file"),
  tar_target(raw_census2022_piped_water,
             "data/raw/census2022/police-district-v1/Police-District/casd-2022-policedist-stata/c2022-pipedwater-policedist-v1.dta",
             format = "file"),
  tar_target(raw_police_districts, "data/raw/boundaries/police_districts/PoliceDistrict.zip",
             format = "file"),
  tar_target(raw_wards_2016, "data/raw/boundaries/mdb_wards_2016/mdb_wards_2016_cra.gpkg",
             format = "file"),
  tar_target(raw_dws_settlements, "data/raw/dws/Umzimvubu_Settlements.gdb", format = "file"),

  # ==========================================================================
  # 20 Prepare
  # ==========================================================================

  # 201 NEX-GDDP hazards, tracked through the ensemble manifest
  script_target("script_201", "R/20_prepare/201_prepare_hazards_nex_gddp.R"),
  tar_target(
    hazards_nex,
    run_script(script_201, "data/processed/hazard/manifests/nex_ensemble_models_used.csv",
               raw_nex_gddp),
    format = "file"
  ),

  # 202 Smallholder presence (SANLC class 41, 20 m)
  script_target("script_202", "R/20_prepare/202_prepare_vulnerability_smallholder_presence.R"),
  tar_target(
    smallholder_presence,
    run_script(script_202,
               file.path(VUL, "smallholder_presence/smallholder_presence_SANLC2022_20m.tif"),
               sanlc_cra),
    format = "file"
  ),

  # 202 alternative: class 41 cropland per household (5 km neighbourhood)
  tar_target(
    rainfed_crop_reliance,
    prepare_cropland_per_household(
      mapdata_zip    = raw_census2011_sal_zip,
      households_dta = raw_census2011_water_source,
      sanlc_file     = sanlc_cra,
      roi            = roi,
      template       = template_grid,
      out_tif        = file.path(VUL, "subsistence_cropland_per_household_1km.tif"),
      out_gpkg_wards = file.path(VUL, "subsistence_cropland_per_household_wards2011.gpkg"),
      out_tif_wards  = file.path(VUL, "subsistence_cropland_per_household_wards2011_1km.tif"),
      crop_classes   = 41,
      radius_m       = 5000,
      min_hh_nbhd    = 50
    ),
    format = "file"
  ),

  # 203 Reliance on agriculture (Census 2011 H14 by 2016 ward)
  tar_target(
    agriculture_reliance,
    prepare_agriculture_reliance(
      union_file     = raw_census2011_activity_union,
      wards_gpkg     = raw_wards_2016,
      sanlc_file     = sanlc_cra,
      template       = template_grid,
      out_gpkg       = file.path(VUL, "agriculture_reliance_census2011_wards2016.gpkg"),
      out_tif        = file.path(VUL, "agriculture_reliance_census2011_1km.tif"),
      out_tif_shares = file.path(VUL, "agriculture_activity_shares_census2011_1km.tif"),
      denominator    = "census2011_sal",
      mapdata_zip    = raw_census2011_sal_zip,
      households_dta = raw_census2011_water_source,
      lm_check_zip   = raw_census2011_agric_sample
    ),
    format = "file"
  ),

  # 204 Grazing capacity
  script_target("script_204", "R/20_prepare/204_prepare_vulnerability_grazing_capacity.R"),
  tar_target(
    grazing_capacity,
    run_script(script_204, file.path(VUL, "grazing_capacity_2016_cra.tif"), raw_grazing),
    format = "file"
  ),

  # 205 Water sources: Census 2022 police districts (rainfall-dependent =
  # primary, unimproved = sensitivity), Census 2011 SALs (unimproved)
  tar_target(
    unimproved_water_2022,
    prepare_unimproved_water_2022(
      source_dta    = raw_census2022_water_source,
      piped_dta     = raw_census2022_piped_water,
      police_zip    = raw_police_districts,
      sanlc_file    = sanlc_cra,
      template      = template_grid,
      out_gpkg      = file.path(VUL, "unimproved_water_police2022.gpkg"),
      out_tif       = file.path(VUL, "unimproved_water_share_census2022_1km.tif"),
      out_tif_rain  = file.path(VUL, "rainfall_dependent_water_share_census2022_1km.tif"),
      out_tif_piped = file.path(VUL, "below_basic_piped_water_census2022_1km.tif")
    ),
    format = "file"
  ),
  tar_target(
    unimproved_water_2011,
    prepare_unimproved_water(
      mapdata_zip   = raw_census2011_sal_zip,
      source_dta    = raw_census2011_water_source,
      piped_dta     = raw_census2011_piped_water,
      sanlc_file    = sanlc_cra,
      roi           = roi,
      template      = template_grid,
      out_gpkg      = file.path(VUL, "unimproved_water_sal2011.gpkg"),
      out_tif       = file.path(VUL, "unimproved_water_share_census2011_1km.tif"),
      out_tif_piped = file.path(VUL, "below_basic_piped_water_census2011_1km.tif"),
      min_hh        = 30
    ),
    format = "file"
  ),

  # 206 Livestock density (GLW4 2020) and stocking pressure
  tar_target(
    livestock_density,
    prepare_livestock_density(
      glw_files    = raw_glw4,
      roi          = roi,
      template     = template_grid,
      grazing      = grazing_capacity,
      out_tlu      = file.path(VUL, "livestock/glw4_2020_tlu_per_km2.tif"),
      out_species  = file.path(VUL, "livestock/glw4_2020_heads_per_km2_species.tif"),
      out_stocking = file.path(VUL, "livestock/stocking_pressure_glw4_2020_daff2016_1km.tif")
    ),
    format = "file"
  ),

  # 207 Old fields: land in lower ecological condition
  tar_target(
    old_fields,
    prepare_old_fields(
      sanlc_file  = sanlc_cra,
      template    = template_grid,
      out_share   = file.path(VUL, "old_fields/old_fields_share_SANLC2022_1km.tif"),
      out_classes = file.path(VUL, "old_fields/old_fields_classes_SANLC2022_20m.tif")
    ),
    format = "file"
  ),

  # 208 Biomass availability (grassland GPP 2022-2024)
  tar_target(
    biomass_gpp,
    prepare_biomass_gpp(raw_gpp, file.path(VUL, "biomass_gpp/grassland_gpp_mean_2022_2024_30m.tif")),
    format = "file"
  ),

  # 209 Alien tree cover (MapWAPS 2025 10 m raster -> 1 km shares)
  tar_target(
    alien_trees,
    prepare_alien_trees(
      mapwaps_tif  = raw_mapwaps[grepl("\\.tif$", raw_mapwaps)],
      lookup       = raw_mapwaps[grepl("lookup", raw_mapwaps)],
      template     = template_grid,
      sanlc_file   = sanlc_cra,
      rivers       = raw_rivers,
      out_cover    = file.path(VUL, "alien_trees/alien_tree_cover_MapWAPS2023_1km.tif"),
      out_invasive = file.path(VUL, "alien_trees/alien_invasive_only_MapWAPS2023_1km.tif"),
      out_riparian = file.path(VUL, "alien_trees/alien_riparian_MapWAPS2023_1km.tif"),
      out_taxa     = file.path(VUL, "alien_trees/alien_taxa_share_MapWAPS2023_1km.tif"),
      out_coverage = file.path(VUL, "alien_trees/mapwaps_coverage_1km.tif")
    ),
    format = "file"
  ),

  # 210 Lower-income areas (Census 2011 income by SAL)
  tar_target(
    low_income_2011,
    prepare_low_income_2011(
      mapdata_zip = raw_census2011_sal_zip,
      income_dta  = raw_census2011_income,
      sanlc_file  = sanlc_cra,
      template    = template_grid,
      out_share   = file.path(VUL, "lower_income/low_income_share_census2011_1km.tif"),
      out_gpkg    = file.path(VUL, "lower_income/low_income_sal2011.gpkg")
    ),
    format = "file"
  ),

  # 211 Settlements: DWS rural households per km2
  tar_target(
    settlement_households,
    prepare_settlement_households(
      dws_gdb     = raw_dws_settlements,
      sanlc_file  = sanlc_cra,
      template    = template_grid,
      out_tif     = file.path(EXP, "settlements/dws_rural_households_per_km2_1km.tif"),
      out_tif_all = file.path(EXP, "settlements/dws_all_households_per_km2_1km.tif"),
      out_gpkg    = file.path(EXP, "settlements/dws_settlements_cra.gpkg"),
      categories  = "Rural"
    ),
    format = "file"
  ),

  # 212 Strategic Water Source Areas, and SWSA degradation for NbS targeting
  tar_target(
    swsa,
    prepare_swsa(
      gpkg         = raw_swsa,
      template     = template_grid,
      out_combined = file.path(EXP, "swsa/swsa_combined_1km.tif"),
      out_surface  = file.path(EXP, "swsa/swsa_surface_1km.tif"),
      out_ground   = file.path(EXP, "swsa/swsa_ground_1km.tif"),
      w_sw         = 0.5
    ),
    format = "file"
  ),
  tar_target(
    swsa_degradation,
    prepare_swsa_degradation(
      surface     = swsa[grepl("swsa_surface", swsa)],
      ground      = swsa[grepl("swsa_ground", swsa)],
      alien       = alien_trees[grepl("alien_tree_cover", alien_trees)],
      gw_pressure = old_fields[grepl("old_fields_share", old_fields)],
      template    = template_grid,
      out_surface = "data/processed/nbs/swsa_degradation/swsa_surface_alien_pressure_1km.tif",
      out_ground  = "data/processed/nbs/swsa_degradation/swsa_ground_oldfields_pressure_1km.tif"
    ),
    format = "file"
  ),

  # ==========================================================================
  # 30 Indicators: one target per active register row (indicator_<id>),
  # written to data/processed/indicators/<component>/<id>_value|score.tif
  # ==========================================================================
  indicator_targets,
  tar_combine(indicator_files, indicator_targets, command = c(!!!.x), format = "file"),
  tar_target(
    indicator_summary,
    summarise_indicators(indicator_files, template_grid,
                         out = "outputs/tables/indicator_summary.csv"),
    format = "file"
  ),

  # ==========================================================================
  # 70 Visualisation
  # ==========================================================================
  script_target("script_701", "R/70_visualisation/701_visualisation_leaflet_map.R"),
  tar_target(framework_xlsx, "../CRA_data_sources_Claude.xlsx", format = "file"),
  tar_target(indicator_register_csv, INDICATOR_REGISTER, format = "file"),
  tar_target(
    indicator_map,
    run_script(script_701, "outputs/maps/umzimvubu_cra_indicators_leaflet.html",
               indicator_files, indicator_register_csv, framework_xlsx,
               catchment_shp, district_shp),
    format = "file"
  ),

  # ==========================================================================
  # 80 Export: map and indicator catalogue published to docs/ for GitHub Pages
  # ==========================================================================
  tar_target(pages_map, export_pages(indicator_map[1]), format = "file"),
  tar_target(
    pages_catalogue,
    export_catalogue(framework_xlsx, indicator_register_csv, indicator_summary),
    format = "file"
  )
)
