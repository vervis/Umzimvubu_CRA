# ============================================================================
# 211_prepare_exposure_settlements.R
#
# Exposure 211 "Settlements": rural households per km2 from the DWS Water
# Services Knowledge System (WSKS) settlement layer (FinYear April 2026,
# Census 2011 modelled and updated annually).
#
# Structure
#   1. Paths and input checks.
#   2. prepare_settlement_households() (R/functions/settlement_functions.R):
#      Rural settlements only; each settlement's TotalHouseholds placed on
#      SANLC 2022 settlement pixels (classes 47-60) inside its polygon at
#      100 m (evenly if none, centroid if tiny), on a grid padded 10 km;
#      summed to the 1 km template and divided by cell area. Template
#      cells without settlements = 0. Writes the indicator, an all-households
#      (rural + urban) sensitivity layer and the settlement table.
#   3. QA: totals, households by LM, distribution, and a preview of the
#      stage-30 score (register row settlement_households: log1p, min-max,
#      trim 0.01).
#
# Not used: DWS population, growth, PopDensity (persons/household) and the
# water-source fields (one profile per LM).
#
# Inputs:  data/raw/dws/Umzimvubu_Settlements.gdb (layer Umzimv_settlements),
#          SANLC 2022 extract, template grid
# Outputs: data/processed/exposure/settlements/
# ============================================================================

library(sf)
library(terra)

source("R/functions/crop_reliance_functions.R")
source("R/functions/water_access_functions.R")
source("R/functions/settlement_functions.R")

gdb   <- "data/raw/dws/Umzimvubu_Settlements.gdb"
sanlc <- "data/processed/landcover/sanlc2022/sanlc2022_umzimvubu_alfred_nzo_extent.tif"
tmpl  <- "data/processed/template/cra_template_1km.tif"
od    <- "data/processed/exposure/settlements"

stopifnot(dir.exists(gdb), file.exists(sanlc), file.exists(tmpl))

out <- prepare_settlement_households(
  dws_gdb     = gdb,
  sanlc_file  = sanlc,
  template    = tmpl,
  out_tif     = file.path(od, "dws_rural_households_per_km2_1km.tif"),
  out_tif_all = file.path(od, "dws_all_households_per_km2_1km.tif"),
  out_gpkg    = file.path(od, "dws_settlements_cra.gpkg"),
  layer       = "Umzimv_settlements",
  categories  = "Rural",
  use_sanlc   = TRUE
)

ru <- rast(out[1])
al <- rast(out[2])
cat("Households in the CRA extent (rural, all):",
    round(unlist(global(c(ru, al), "sum", na.rm = TRUE))), "\n")

s  <- st_drop_geometry(st_read(out[3], quiet = TRUE))
lm <- aggregate(TotalHouseholds ~ LMName + Category, s, sum)
print(lm[order(lm$Category, -lm$TotalHouseholds), ], row.names = FALSE)

v <- values(ru, na.rm = TRUE)
cat("Template cells:", length(v), "| with rural households:", sum(v > 0),
    sprintf("(%.0f%%)", 100 * mean(v > 0)), "\n")
print(round(quantile(v[v > 0], c(0.1, 0.5, 0.9, 0.99, 1)), 1))

lx <- log1p(ru)
hi <- quantile(values(lx, na.rm = TRUE), 0.99)
sc <- clamp(lx / hi, 0, 1)
cat("Score among populated cells (p10, p50, p90):",
    round(quantile(values(sc)[values(ru) > 0], c(0.1, 0.5, 0.9), na.rm = TRUE), 2),
    "| threshold for score 1:", round(expm1(hi)), "households/km2\n")

par(mfrow = c(1, 2))
plot(lx, main = "log1p rural households per km2")
plot(sc, main = "Score preview (log1p, min-max, trim 0.01)")
