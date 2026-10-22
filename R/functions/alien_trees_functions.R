# ============================================================================
# R/functions/alien_trees_functions.R
#
# Vulnerability 209 "Alien tree cover": functions used by
# 110_download_vulnerability_mapwaps.R, 209_prepare_vulnerability_alien_trees.R
# and _targets.R (targets raw_mapwaps and alien_trees).
#
# Source: MapWAPS Invasive Alien Plant Map for the uMzimvubu Catchment
# (Skosana, Cogill & Rebelo, Stellenbosch University; SUNScholarData,
# doi:10.25413/sun.25050401). Random Forest classification of Sentinel-2
# (15 May 2023). Used: the 19 Aug 2025 release, a 10 m class raster
# (MAPWAPS_uMzimvubu_19082025.tif, WGS 84 / UTM 35S, 18 classes, 0 =
# background). Figshare replaced it in Aug 2026 with a 2 GB shapefile of 15
# class-level multipolygons that GDAL cannot read in practical time, and keeps
# no older files, so the tif is held locally and not downloaded (archived
# shapefile release: data/raw/mapwaps/Newer_shp_version). Covers the
# catchment only, so Alfred Nzo cells outside it are NA. Plantations are
# included in the main layer: pine and gum use the same water whether planted
# or invading.
#
# Constants
#   MAPWAPS_DOI        Figshare record
#   MAPWAPS_CLASSES    2025 legend (from the authors' QGIS style); alien = 1,
#                      plantation_taxon = 1 for pine and gum
#   SANLC_PLANTATION_CLASSES  SANLC 2022 planted forest (5 dense, 6 open,
#                      7 clear-felled)
#
# Functions
#   register_mapwaps()    checks the local tif, writes the lookup, a manifest
#                         (bytes, md5) and a .url; returns the tif path
#   write_mapwaps_lookup(), read_mapwaps_lookup()
#   prepare_alien_trees() 10 m classes -> 100 m counts -> 1 km shares on the
#                         template; optional invasive-only and riparian layers
# ============================================================================


MAPWAPS_DOI     <- "https://doi.org/10.25413/sun.25050401"
MAPWAPS_TIF     <- "MAPWAPS_uMzimvubu_19082025.tif"
MAPWAPS_RELEASE <- "2025-08-19"

MAPWAPS_CLASSES <- data.frame(
  value = 0:18,
  label = c("Background", "Black/Green Wattle", "Gum",
            "Other Invasive Alien Plants", "Pine", "Poplar", "Silver Wattle",
            "Bare Ground", "Dryland Agriculture", "Grassland",
            "Indigenous Bush (Other)", "Indigenous Bush (Ouhoudt)",
            "Indigenous Bush (Vachellia)", "Irrigated Agriculture", "Maize",
            "Burnt", "Urban", "Water", "Wetland"),
  taxon = c(NA, "wattle_black", "gum", "other_alien", "pine", "poplar",
            "wattle_silver", rep(NA, 12)),
  alien = c(0, rep(1, 6), rep(0, 12)),
  plantation_taxon = c(0, 0, 1, 0, 1, 0, 0, rep(0, 12))
)

SANLC_PLANTATION_CLASSES <- 5:7


register_mapwaps <- function(dir = "data/raw/mapwaps", tif = MAPWAPS_TIF) {
  f <- file.path(dir, tif)
  if (!file.exists(f)) {
    stop(f, " not found. The 2025 MapWAPS raster is no longer on Figshare: ",
         "restore it from Dropbox (Deleted files / version history) or the authors.")
  }
  write_mapwaps_lookup(file.path(dir, "mapwaps_class_lookup.csv"))
  writeLines(c("[InternetShortcut]", paste0("URL=", MAPWAPS_DOI)),
             file.path(dir, "mapwaps_umzimvubu.url"))
  utils::write.csv(data.frame(
    file = tif, bytes = file.size(f), md5 = unname(tools::md5sum(f)),
    release = MAPWAPS_RELEASE, doi = MAPWAPS_DOI, licence = "CC BY 4.0",
    note = "Superseded on Figshare (Aug 2026, shapefile only); held locally",
    checked = format(Sys.Date())
  ), file.path(dir, "mapwaps_manifest.csv"), row.names = FALSE)
  f
}


write_mapwaps_lookup <- function(out, classes = MAPWAPS_CLASSES) {
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(classes, out, row.names = FALSE, na = "")
  out
}

read_mapwaps_lookup <- function(lookup, r = NULL) {
  lk <- utils::read.csv(lookup, stringsAsFactors = FALSE)
  if (!is.null(r)) {
    mm <- terra::minmax(r, compute = TRUE)
    if (any(is.na(mm)) || mm[1] < min(lk$value) || mm[2] > max(lk$value)) {
      stop("MapWAPS values (", paste(mm, collapse = "-"), ") fall outside the ",
           "lookup (", min(lk$value), "-", max(lk$value), "): update ", basename(lookup), ".")
    }
  }
  lk
}


# ----------------------------------------------------------------------------
# prepare_alien_trees(): MapWAPS classes -> alien tree shares on 1 km grid
# ----------------------------------------------------------------------------
#   mapwaps_tif : 10 m MapWAPS class raster (2025 release)
#   lookup      : class lookup CSV (value, label, taxon, alien, plantation_taxon)
#   template    : 1 km CRA template
#   sanlc_file  : SANLC 2022 (optional). Pine/gum pixels on SANLC planted
#                 forest are dropped for the invasive-only companion.
#   rivers      : river lines (optional, e.g. 111 NFEPA rivers) for the
#                 riparian companion: alien pixels within riparian_m of a
#                 river, as a share of the mapped cell (0 = no river).
#                 riparian_m = 100 allows for the ~1:500 000 line accuracy.
#   min_cov     : minimum share of a 1 km cell mapped by MapWAPS; below it
#                 the cell is NA (catchment edge).
# Counts are summed to 100 m in UTM (exact), then averaged onto the template,
# so share = alien pixels / mapped pixels.

prepare_alien_trees <- function(mapwaps_tif, lookup, template,
                                out_cover, out_taxa, out_coverage,
                                sanlc_file = NULL, out_invasive = NULL,
                                rivers = NULL, out_riparian = NULL,
                                plantation_classes = SANLC_PLANTATION_CLASSES,
                                riparian_m = 100, min_cov = 0.5, fact = 10) {
  m    <- terra::rast(mapwaps_tif)[[1]]
  lk   <- read_mapwaps_lookup(lookup, m)
  tmpl <- terra::rast(template)
  al   <- lk[lk$alien == 1, ]

  rc <- function(v) terra::classify(m, cbind(lk$value, v), others = 0)
  lay <- list(
    mapped = rc(as.numeric(lk$value > 0)),
    alien  = rc(lk$alien)
  )
  for (i in seq_len(nrow(al))) {
    lay[[al$taxon[i]]] <- rc(as.numeric(lk$value == al$value[i]))
  }

  if (!is.null(sanlc_file)) {
    s  <- terra::rast(sanlc_file)[[1]]
    pl <- terra::ifel(s %in% plantation_classes, 1, 0)
    pl <- terra::project(pl, m, method = "near")
    pt <- rc(lk$plantation_taxon)
    lay$alien_invasive <- lay$alien - pt * terra::subst(pl, NA, 0)
  }

  if (!is.null(rivers)) {
    e  <- as.vector(terra::ext(m))
    bb <- sf::st_bbox(c(xmin = e[[1]], xmax = e[[2]], ymin = e[[3]], ymax = e[[4]]),
                      crs = sf::st_crs(terra::crs(m)))
    rv <- sf::st_read(rivers, quiet = TRUE) |>
      sf::st_zm() |>
      sf::st_transform(sf::st_crs(bb)) |>
      sf::st_crop(bb)
    rb <- sf::st_buffer(sf::st_union(rv), riparian_m)
    rz <- terra::rasterize(terra::vect(rb), m, field = 1, background = 0)
    lay$alien_riparian <- lay$alien * rz
  }

  cnt <- terra::aggregate(terra::rast(lay), fact = fact, fun = "sum", na.rm = TRUE)
  km  <- terra::project(cnt, tmpl, method = "average")
  km  <- terra::mask(km, tmpl)

  px  <- fact^2
  cov <- km[["mapped"]] / px
  ok  <- cov >= min_cov
  shr <- function(x) terra::mask(x / km[["mapped"]], ok, maskvalues = c(0, NA))

  cover <- shr(km[["alien"]])
  names(cover) <- "alien_tree_cover"
  taxa <- shr(km[[al$taxon]])
  names(taxa) <- al$taxon
  names(cov) <- "mapwaps_coverage"

  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  wr <- function(x, f) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    terra::writeRaster(x, f, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
    f
  }
  out <- c(cover = wr(cover, out_cover))
  if (!is.null(sanlc_file) && !is.null(out_invasive)) {
    inv <- shr(km[["alien_invasive"]])
    names(inv) <- "alien_invasive_only"
    out <- c(out, invasive = wr(inv, out_invasive))
  }
  if (!is.null(rivers) && !is.null(out_riparian)) {
    rip <- shr(km[["alien_riparian"]])
    names(rip) <- "alien_riparian"
    out <- c(out, riparian = wr(rip, out_riparian))
  }
  c(out, taxa = wr(taxa, out_taxa), coverage = wr(cov, out_coverage))
}
