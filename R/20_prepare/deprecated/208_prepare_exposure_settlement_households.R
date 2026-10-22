# ============================================================================
# 208_prepare_exposure_settlement_households.R
#
# Exposure 208 "Settlements": rural households per km2.
#
# Exposure here is the presence of the communities whose water and farming
# could be hit. The layer counts households, not settlement area: SANLC
# settlement area and DWS households agree on where people live, but the
# number of households per ha of settlement varies about 15-fold between
# cells, and SANLC cannot separate villages from towns.
#
# Source: DWS Water Services Knowledge System (WSKS) settlement layer,
# FinYear April 2026 ("Census 2011 modelled and updated", annual). One
# polygon per settlement with a household count and a Rural/Urban category.
# Received from DWS in Oct 2026: data/raw/dws/Umzimvubu_Settlements.gdb.
#
# Method (prepare_settlement_households() in R/functions/settlement_functions.R)
#   1. Rural settlements only (Category == "Rural"); the target communities
#      are rural smallholders. An all-households companion layer is written
#      for a sensitivity run.
#   2. Each settlement's households are placed on the SANLC 2022 settlement
#      pixels (classes 47-60) inside its polygon at 100 m, evenly if it has
#      none, then summed to the 1 km template grid and divided by cell area.
#   3. Study-area cells without a DWS settlement are 0, not NA.
#
# Stage 30 (register row settlement_households): log1p transform, then
# min-max with trim 0.01, so the many zeros stay 0 and a few dense villages
# do not squash everything else. Test run 9 Oct 2026: 40% of study-area cells
# have households; median 17, p90 84, p99 216, max ~1,230 households/km2;
# ~300,000 rural households in the study area.
#
# Known gaps: DWS misses ~10% of Census 2022 households, mostly in towns and
# on commercial farms (Ngangelizwe, Mthatha, Maloti, Ugie, Richmond-KZN).
# Do not use the DWS population, growth or "PopDensity" fields.
#
# Inputs
#   data/raw/dws/Umzimvubu_Settlements.gdb (layer Umzimv_settlements)
#   data/processed/landcover/sanlc2022/sanlc2022_umzimvubu_alfred_nzo_extent.tif
#   data/processed/template/cra_template_1km.tif (run
#     targets::tar_make(template_grid) first if it does not exist)
#
# Outputs (data/processed/exposure/settlements/)
#   dws_rural_households_per_km2_1km.tif  the indicator
#   dws_all_households_per_km2_1km.tif    rural + urban (sensitivity)
#   dws_settlements_cra.gpkg              settlements touching the study area
#
# The pipeline runs the same function through _targets.R (target
# settlement_households). Needs the sf and terra packages.
# ============================================================================

library(sf)
library(terra)

source("R/functions/crop_reliance_functions.R")    # spread_households()
source("R/functions/water_access_functions.R")     # SANLC_SETTLEMENT_CLASSES
source("R/functions/settlement_functions.R")

# 1. Paths --------------------------------------------------------------------

dws_gdb    <- "data/raw/dws/Umzimvubu_Settlements.gdb"
sanlc_file <- "data/processed/landcover/sanlc2022/sanlc2022_umzimvubu_alfred_nzo_extent.tif"
template   <- "data/processed/template/cra_template_1km.tif"
out_dir    <- "data/processed/exposure/settlements"

stopifnot(dir.exists(dws_gdb), file.exists(sanlc_file), file.exists(template))


# 2. Rural households per km2 -------------------------------------------------

out <- prepare_settlement_households(
  dws_gdb     = dws_gdb,
  sanlc_file  = sanlc_file,
  template    = template,
  out_tif     = file.path(out_dir, "dws_rural_households_per_km2_1km.tif"),
  out_tif_all = file.path(out_dir, "dws_all_households_per_km2_1km.tif"),
  out_gpkg    = file.path(out_dir, "dws_settlements_cra.gpkg"),
  categories  = "Rural",   # c("Rural", "Urban") for all households
  use_sanlc   = TRUE       # FALSE: spread evenly over each polygon
)


# 3. QA -----------------------------------------------------------------------

rural <- rast(out[1])
allhh <- rast(out[2])

# Households in the study area (cells are 1 km2, so the sum is a count)
print(global(c(rural, allhh), "sum", na.rm = TRUE))

# Rural households by local municipality (settlements touching the study area)
set <- st_read(out[3], quiet = TRUE)
lm <- aggregate(TotalHouseholds ~ LMName + Category, data = st_drop_geometry(set), sum)
print(lm[order(lm$Category, -lm$TotalHouseholds), ])

# Preview of the score stage 30 will make (log1p, min-max, trim 0.01)
lx <- log1p(rural)
hi <- global(lx, function(x) quantile(x, 0.99, na.rm = TRUE))[[1]]
plot(clamp(lx / hi, 0, 1), main = "Rural households per km2 (score preview)")
