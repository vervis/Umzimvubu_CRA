# 212 Strategic Water Source Areas: tier-weighted surface (2021) + groundwater (2017) shares, combined exposure, SWSA degradation (NbS targeting)
# Needs 207 (old fields) and 209 (alien trees) outputs for the degradation layers.

library(terra)

source("R/functions/swsa_functions.R")

gpkg     <- "data/raw/swsa/Strategic_Water_Source_Areas_CRA.gpkg"
template <- "data/processed/template/cra_template_1km.tif"
alien    <- "data/processed/vulnerability/alien_trees/alien_tree_cover_MapWAPS2023_1km.tif"
gw_pres  <- "data/processed/vulnerability/old_fields/old_fields_share_SANLC2022_1km.tif"
out_dir  <- "data/processed/exposure/swsa"
deg_dir  <- "data/processed/nbs/swsa_degradation"

stopifnot(file.exists(gpkg), file.exists(template))

out <- prepare_swsa(
  gpkg         = gpkg,
  template     = template,
  out_combined = file.path(out_dir, "swsa_combined_1km.tif"),
  out_surface  = file.path(out_dir, "swsa_surface_1km.tif"),
  out_ground   = file.path(out_dir, "swsa_ground_1km.tif")
)

tmpl <- rast(template)
x <- cbind(e  = values(rast(out[["combined"]]), mat = FALSE),
           sw = values(rast(out[["surface"]]),  mat = FALSE),
           gw = values(rast(out[["ground"]]),   mat = FALSE))
x <- x[!is.na(values(tmpl, mat = FALSE)), ]
cat("Study-area cells:", nrow(x), "\n")
print(round(colMeans(x > 0), 3))
print(round(apply(x, 2, quantile, c(0.25, 0.5, 0.75, 0.9)), 2))
print(table(cut(x[, "e"], c(-Inf, 0, 0.25, 0.5, 0.75, 1), right = TRUE)))

if (file.exists(alien) && file.exists(gw_pres)) {
  deg <- prepare_swsa_degradation(
    surface     = out[["surface"]],
    ground      = out[["ground"]],
    alien       = alien,
    gw_pressure = gw_pres,
    template    = template,
    out_surface = file.path(deg_dir, "swsa_surface_alien_pressure_1km.tif"),
    out_ground  = file.path(deg_dir, "swsa_ground_oldfields_pressure_1km.tif")
  )
  for (p in list(c("surface", alien), c("ground", gw_pres))) {
    s <- values(rast(out[[p[1]]]), mat = FALSE)
    v <- values(rast(p[2]), mat = FALSE)
    ok <- !is.na(values(tmpl, mat = FALSE)) & !is.na(v)
    cat(p[1], "| mean pressure in national SWSA:", round(mean(v[ok & s > 0.5]), 3),
        "| outside:", round(mean(v[ok & s == 0]), 3),
        "| share of pressure in national SWSA:", round(sum(v[ok & s > 0.5]) / sum(v[ok]), 2),
        "| Spearman:", round(cor(s[ok], v[ok], method = "spearman"), 2), "\n")
  }
  par(mfrow = c(2, 2))
  plot(rast(out[["surface"]]), main = "Surface SWSA (tier-weighted share)")
  plot(rast(out[["ground"]]), main = "Groundwater SWSA (tier-weighted share)")
  plot(rast(deg[["surface"]]), main = "Surface SWSA x alien tree cover")
  plot(rast(deg[["ground"]]), main = "Groundwater SWSA x old fields")
} else {
  message("Degradation layers skipped: run 207 and 209 first.")
  par(mfrow = c(1, 3))
  plot(rast(out[["surface"]]), main = "Surface SWSA")
  plot(rast(out[["ground"]]), main = "Groundwater SWSA")
  plot(rast(out[["combined"]]), main = "212 combined")
}
