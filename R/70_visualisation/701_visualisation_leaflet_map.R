# ============================================================================
# 701_visualisation_leaflet_map.R
#
# Interactive Leaflet map of the CRA layers: stage-30 indicators
# (data/processed/indicators/<component>/<id>_value|score.tif), organised by
# the framework hierarchy in CRA_data_sources_Claude.xlsx (sheet
# Indicators: Sub-risk > Component > Factor > Indicator), and the stage
# 40-50 composite indices (hazard, vulnerability, exposure, risk) with
# their GCM ensemble statistics.
#   1. Settings, boundaries and CRA extent.
#   2. Framework: xlsx hierarchy joined to the register on Register ID; rows
#      without stage-30 outputs are listed but disabled (register fallback).
#   3. Climate model data on the native NEX grid (201): per-GCM values of
#      every NEX indicator and period, per-GCM significance and the IPCC AR6
#      robustness class of the change (1 robust, 2 no robust signal,
#      3 conflicting). Used client-side for the change view (diverging
#      colours by change towards higher risk, hatching: none = robust,
#      / = no change or no robust signal, X = conflicting) and for the
#      location query.
#   4. Location query data: static indicator values on the 1 km template
#      (Float32, base64), normalisation bounds (stage 30 _bounds.csv),
#      component weights and risk settings (stages 40-50). A click
#      reproduces the pipeline at that 1 km cell: every GCM's hazard
#      indicators are interpolated bilinearly from the native grid (as
#      terra::project does), scored with the shared bounds, combined into
#      hazard and risk per GCM, and summarised (median, P10-P90, share of
#      GCMs at or above the high-risk threshold).
#   5. Rasters (PNG overlays): indicator values/scores per period (viridis
#      pooled across periods; scores fixed 0-1 YlOrRd) and composite
#      indices (median, P10, P90, % of GCMs high risk).
#   6. Map + panel (onRender JS): layer tree with composite indices,
#      Values / Scores / Change from historical, separate scenario (SSP) and
#      period selectors, ensemble statistic for composites, hatching toggle,
#      opacity, boundaries, basemaps, legend, About & data sources, and a
#      location query tool (toggle button) with a results panel and a
#      per-GCM chart (historical + every SSP for the selected period).
#      Opens on the layer named in the URL hash (index.html#<id>).
#      Saved to outputs/maps/umzimvubu_cra_indicators_leaflet.html.
# ============================================================================

library(sf)
library(terra)
library(leaflet)
library(htmlwidgets)
library(htmltools)
library(readxl)
library(jsonlite)

source("R/functions/pipeline_helpers.R")
source("R/functions/indicator_functions.R")

reg_csv  <- "R/30_indicators/indicator_register.csv"
fw_xlsx  <- "../CRA_data_sources_Claude.xlsx"
ind_dir  <- "data/processed/indicators"
cmp_dir  <- "data/processed/components"
rsk_dir  <- "data/processed/risk"
man_csv  <- "data/processed/hazard/manifests/nex_ensemble_models_used.csv"
out_html <- "outputs/maps/umzimvubu_cra_indicators_leaflet.html"
opac     <- 0.8
repo_url <- "https://github.com/vervis/Umzimvubu_CRA"
cat_url  <- "https://vervis.github.io/Umzimvubu_CRA/indicators.html"
about    <- paste(
  "Indicators and composite risk for SR1 of the Umzimvubu Climate Risk Assessment",
  "(ALBATROSS, Umzimvubu Hub), at 1 km resolution. Following the IPCC AR6 risk",
  "concept, risk combines hazard, vulnerability and exposure; scores are scaled",
  "0-1 (1 = higher risk) across all periods, GCMs and scenarios. Climate",
  "hazards come from the NEX-GDDP-CMIP6 ensemble (~25 km): maps show the",
  "ensemble median, the change view shows the median change with IPCC AR6",
  "robustness hatching, and the location query shows the spread across GCMs.",
  "Draft research output: indicative only, not for site-level decisions."
)
ssp_lab <- c(ssp126 = "SSP1-2.6", ssp245 = "SSP2-4.5", ssp370 = "SSP3-7.0",
             ssp585 = "SSP5-8.5")
dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)

rd <- function(f) st_transform(st_zm(st_make_valid(st_read(f, quiet = TRUE))), 4326)
b64 <- function(x, size) base64_enc(writeBin(as.integer(x), raw(), size = size, endian = "little"))

bnd <- list(
  list(g = "Umzimvubu catchment", f = CRA_BOUNDARIES[["catchment"]],
       col = "#1f78b4", wt = 3, dash = "", on = TRUE),
  list(g = "Alfred Nzo District Municipality", f = CRA_BOUNDARIES[["district"]],
       col = "#e66101", wt = 3, dash = "8 6", on = TRUE),
  list(g = "Matatiele Local Municipality", f = "data/raw/boundaries/Matatiele.shp",
       col = "#7b3294", wt = 2, dash = "2 5", on = FALSE)
)
bnd <- Filter(function(b) file.exists(b$f), bnd)

e   <- as.vector(cra_extent())
box <- st_as_sfc(st_bbox(c(xmin = e[[1]], xmax = e[[2]], ymin = e[[3]], ymax = e[[4]]),
                         crs = st_crs(CRA_CRS)))
box <- st_transform(st_segmentize(box, CRA_RES), 4326)
roi <- st_bbox(box)


# 2. Framework -----------------------------------------------------------------

reg <- read.csv(reg_csv, stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (!"weight" %in% names(reg)) reg$weight <- 1

if (file.exists(fw_xlsx)) {
  x  <- as.data.frame(read_excel(fw_xlsx, sheet = "Indicators", skip = 1))
  x  <- x[!is.na(x[["Register ID"]]), ]
  fw <- data.frame(id = x[["Register ID"]], subrisk = x[["Sub-risk"]],
                   component = x[["Component"]], factor = x[["Factor"]],
                   q = x[["Question answered"]], role = x[["Role"]],
                   ds = x[["Dataset"]], pv = x[["Provider"]],
                   yr = x[["Data year / coverage"]], lic = x[["Access / licence"]],
                   url = x[["Data URL"]], cite = x[["Documentation / citation"]])
  lk  <- as.data.frame(read_excel(fw_xlsx, sheet = "Lookups", col_names = FALSE,
                                  .name_repair = "minimal"))
  lk  <- lk[grepl("^SR\\d", lk[[1]]), ]
  srl <- setNames(as.character(lk[[2]]), lk[[1]])
} else {
  fw <- data.frame(id = reg$id, subrisk = reg$subrisk,
                   component = tools::toTitleCase(reg$component),
                   factor = reg$label, q = NA, role = NA, ds = reg$source,
                   pv = NA, yr = NA, lic = NA, url = NA, cite = NA)
  srl <- character()
}

k  <- match(fw$id, reg$id)
fw <- fw[!is.na(k), ]
k  <- k[!is.na(k)]
ind <- cbind(fw, reg[k, c("label", "units", "temporal", "direction", "src")],
             cdir = tolower(reg$component[k]))
ind$label <- ifelse(is.na(ind$label), ind$id, ind$label)
ind$fv <- file.path(ind_dir, ind$cdir, paste0(ind$id, "_value.tif"))
ind$fs <- file.path(ind_dir, ind$cdir, paste0(ind$id, "_score.tif"))
ind$fb <- file.path(ind_dir, ind$cdir, paste0(ind$id, "_bounds.csv"))
ind$built <- file.exists(ind$fv) & file.exists(ind$fs)
if (!any(ind$built)) stop("No indicator outputs found in ", ind_dir)

pinfo <- function(p) {
  if (p == "baseline") return(list(hz = "baseline", ssp = ""))
  s <- regmatches(p, regexpr("ssp\\d{3}", p))
  if (!length(s)) return(list(hz = "hist", ssp = ""))
  list(hz = sub("_.*", "", p), ssp = s)
}
plab <- function(p) {
  if (p == "baseline") return("Baseline (static layer)")
  y <- sub("_", "-", regmatches(p, regexpr("\\d{4}_\\d{4}", p)))
  s <- regmatches(p, regexpr("ssp\\d{3}", p))
  if (!length(s)) return(paste("Historical", y))
  h <- c(near = "Near", mid = "Mid", late = "Late")[sub("_.*", "", p)]
  sprintf("%s %s, %s", h, y, ssp_lab[[s]])
}


# 3. Climate model data (native NEX grid) -------------------------------------

nexi <- ind[ind$built & ind$temporal == "nex", ]
nx <- NULL
if (nrow(nexi) && file.exists(man_csv)) {
  models <- sort(unique(read.csv(man_csv)$model))
  nper <- names(rast(nexi$fv[1]))
  g0 <- rast(gcm_sources(nexi$src[1], nper[1], models)[[1]][1])
  nx <- list(x0 = xmin(g0), y1 = ymax(g0), dx = res(g0)[1], dy = res(g0)[2],
             nc = ncol(g0), nr = nrow(g0), models = models, ind = list())
  for (i in seq_len(nrow(nexi))) {
    gs <- gcm_sources(nexi$src[i], nper, models)
    v <- lapply(gs, function(fs) signif(as.vector(values(rast(fs))), 5))
    sg <- cl <- list()
    for (p in nper[!grepl("^hist", nper)]) {
      fsg <- sub("ensemble_median", "significant", gsub("{period}", p, nexi$src[i], fixed = TRUE))
      fcl <- sub("ensemble_median", "robust_class", gsub("{period}", p, nexi$src[i], fixed = TRUE))
      if (file.exists(fsg)) {
        z <- as.vector(values(rast(fsg)[[models]]))
        sg[[p]] <- paste(ifelse(is.na(z), "-", z), collapse = "")
      }
      if (file.exists(fcl)) {
        z <- values(rast(fcl), mat = FALSE)
        cl[[p]] <- I(ifelse(is.na(z), 0L, as.integer(z)))
      }
    }
    nx$ind[[nexi$id[i]]] <- list(
      rel = grepl("annual_rainfall|pet_", nexi$src[i]), v = v, sig = sg, cls = cl
    )
  }
}


# 4. Location query data (1 km template) ------------------------------------

bnds <- lapply(setNames(ind$fb[ind$built], ind$id[ind$built]), function(f)
  if (file.exists(f)) as.list(read.csv(f)) else NULL)

sti <- ind[ind$built & ind$temporal == "static", ]
t0  <- rast(sti$fv[1])
qd  <- list(x0 = xmin(t0), y1 = ymax(t0), res = res(t0)[1], nc = ncol(t0), nr = nrow(t0),
            st = list())
for (i in seq_len(nrow(sti))) {
  v  <- values(rast(sti$fv[i])[[1]], mat = FALSE)
  b  <- bnds[[sti$id[i]]]
  z  <- list(v = base64_enc(writeBin(as.numeric(ifelse(is.na(v), NaN, v)), raw(), size = 4,
                                     endian = "little")))
  if (is.null(b) || b$normalise != "minmax") {
    s <- values(rast(sti$fs[i])[[1]], mat = FALSE)
    z$s <- b64(ifelse(is.na(s), 255L, round(s * 250)), 1)
  }
  qd$st[[sti$id[i]]] <- z
}
pp <- regmatches(CRA_CRS, gregexpr("\\+(lat_0|lat_1|lat_2|lon_0)=[-0-9.]+", CRA_CRS))[[1]]
aea <- setNames(as.list(as.numeric(sub(".*=", "", pp))), sub("\\+(.*)=.*", "\\1", pp))

cmp <- NULL
cf  <- list.files(cmp_dir, "_components\\.csv$", full.names = TRUE)
rf  <- list.files(rsk_dir, "_risk_settings\\.csv$", full.names = TRUE)
if (length(cf) && length(rf)) {
  ct <- read.csv(cf[1], stringsAsFactors = FALSE)
  rs <- read.csv(rf[1], stringsAsFactors = FALSE)
  sr <- ct$subrisk[1]
  cmp <- list(
    subrisk = sr, method = rs$method, thr = rs$threshold, q_high = rs$q_high,
    min_w = ct$min_w[1],
    comps = lapply(setNames(seq_len(nrow(ct)), ct$component), function(j) list(
      ids = I(strsplit(ct$indicators[j], ";")[[1]]),
      w = I(as.numeric(strsplit(ct$weights[j], ";")[[1]])),
      dyn = ct$dynamic[j], file = ct$file[j]
    )),
    risk = file.path(rsk_dir, paste0(sr, "_risk_ensemble.tif"))
  )
}


# 5. Rasters + 6. Map ------------------------------------------------------------

bm <- c("Topographic (Esri)" = "Esri.WorldTopoMap",
        "Satellite (Esri)" = "Esri.WorldImagery")

m <- leaflet(options = leafletOptions(preferCanvas = TRUE))
for (b in names(bm)) m <- addProviderTiles(m, bm[[b]], group = b)

m <- addPolygons(m, data = box, group = "CRA extent", fill = FALSE,
                 color = "#555555", weight = 1.5, dashArray = "4 4")
for (b in bnd) {
  m <- addPolygons(m, data = rd(b$f), group = b$g, fill = FALSE, color = b$col,
                   weight = b$wt, opacity = 1, label = b$g,
                   dashArray = if (nzchar(b$dash)) b$dash else NULL)
}

lyr <- character()
leg <- list()
per <- list()
add_layers <- function(m, r, key, pal_name, rg) {
  if (!all(is.finite(rg))) return(m)
  if (diff(rg) == 0) rg <- rg + c(-0.5, 0.5)
  pal <- colorNumeric(pal_name, rg, na.color = "transparent")
  for (j in seq_len(nlyr(r))) {
    g <- paste(key, names(r)[j], sep = "|")
    m <- addRasterImage(m, r[[j]], colors = pal, opacity = opac,
                        project = TRUE, method = "ngb", group = g)
    lyr <<- c(lyr, g)
  }
  leg[[key]] <<- list(cols = pal(seq(rg[1], rg[2], length.out = 9)), min = rg[1], max = rg[2])
  m
}
rng <- function(r) range(values(r, mat = FALSE), na.rm = TRUE)

for (i in which(ind$built)) {
  for (ty in c("value", "score")) {
    r <- rast(ind[[if (ty == "value") "fv" else "fs"]][i])
    m <- add_layers(m, r, paste(ind$id[i], ty, sep = "|"),
                    if (ty == "score") "YlOrRd" else "viridis",
                    if (ty == "score") c(0, 1) else rng(r))
    per[[ind$id[i]]] <- names(r)
  }
}
ind$built <- ind$id %in% names(per)

citems <- list()
if (!is.null(cmp)) {
  clab <- c(risk = "Climate risk (H x V x E)", hazard = "Hazard index",
            vulnerability = "Vulnerability index", exposure = "Exposure index")
  slab <- c(median = "Median of GCMs", p10 = "P10 of GCMs", p90 = "P90 of GCMs",
            high_share = "% of GCMs with high risk", index = "Index")
  add_cmp <- function(m, key, r, stats) {
    pp <- unique(sub(paste0("_(", paste(stats, collapse = "|"), ")$"), "", names(r)))
    for (s in stats) {
      rr <- r[[paste0(pp, "_", s)]]
      names(rr) <- pp
      m <- add_layers(m, rr, paste("cmp", key, s, sep = "|"),
                      if (s == "high_share") "Purples" else "YlOrRd",
                      if (s == "high_share") c(0, 100) else
                        range(values(r[[!grepl("_high_share$", names(r))]], mat = FALSE), na.rm = TRUE))
    }
    citems[[length(citems) + 1]] <<- list(key = key, label = clab[[key]], stats = I(stats),
                                          slabs = as.list(slab[stats]), periods = I(pp))
    m
  }
  rr <- rast(cmp$risk)
  m <- add_cmp(m, "risk", rr, c("median", "p10", "p90", "high_share"))
  for (cp in c("hazard", "vulnerability", "exposure")) {
    z <- cmp$comps[[cp]]
    if (is.null(z)) next
    r <- rast(z$file)
    if (z$dyn) {
      m <- add_cmp(m, cp, r, c("median", "p10", "p90"))
    } else {
      names(r) <- "baseline_index"
      m <- add_cmp(m, cp, r, "index")
    }
  }
}

pers <- unique(c(unlist(per), unlist(lapply(citems, `[[`, "periods"))))
pers <- c(setdiff(pers, "baseline"), intersect(pers, "baseline"))
pk <- list()
for (p in pers) {
  z <- pinfo(p)
  if (z$hz == "hist") pk$hist <- p
  if (nzchar(z$ssp)) pk[[paste(z$hz, z$ssp, sep = "|")]] <- p
}
ssps <- intersect(names(ssp_lab), unique(unlist(lapply(pers, function(p) pinfo(p)$ssp))))
hzs  <- intersect(c("near", "mid", "late"), unique(unlist(lapply(pers, function(p) pinfo(p)$hz))))
hzl  <- vapply(hzs, function(h) {
  p <- pk[[paste(h, ssps[1], sep = "|")]]
  paste(tools::toTitleCase(h), sub("_", "-", regmatches(p, regexpr("\\d{4}_\\d{4}", p))))
}, "")
hl <- if (!is.null(pk$hist)) paste("Historical", sub("_", "-", regmatches(pk$hist, regexpr("\\d{4}_\\d{4}", pk$hist)))) else ""
dflt <- ind[ind$built, ][1, ]

urls <- function(u) {
  if (is.na(u)) return(character())
  u <- trimws(unlist(strsplit(u, "[;[:space:]]+")))
  u[grepl("^https?://", u)]
}

cmp_ord <- c("Hazard", "Vulnerability", "Exposure")
tree <- lapply(unique(ind$subrisk), function(s) {
  a  <- ind[ind$subrisk == s, ]
  cs <- unique(a$component)
  cs <- c(intersect(cmp_ord, cs), setdiff(cs, cmp_ord))
  list(code = s, desc = if (s %in% names(srl)) srl[[s]] else "",
       comps = lapply(cs, function(cc) {
         b <- a[a$component == cc, ]
         list(name = cc, factors = lapply(unique(b$factor), function(fc) {
           d <- b[b$factor == fc, ]
           list(name = fc, inds = lapply(seq_len(nrow(d)), function(j) list(
             id = d$id[j], label = d$label[j], q = d$q[j], role = d$role[j],
             units = d$units[j], dir = d$direction[j], built = d$built[j],
             comp = d$cdir[j], nex = d$temporal[j] == "nex",
             periods = I(if (d$built[j]) per[[d$id[j]]] else character()),
             src = list(ds = d$ds[j], pv = d$pv[j], yr = d$yr[j], lic = d$lic[j],
                        url = I(urls(d$url[j])), cite = d$cite[j])
           )))
         }))
       }))
})

cfg <- list(
  title = "Umzimvubu CRA indicators", tree = tree, layers = I(lyr), legends = leg,
  periods = I(pers), plabs = setNames(lapply(pers, plab), pers), pk = pk,
  ssps = I(ssps), ssplabs = as.list(ssp_lab[ssps]), hzs = I(hzs), hzlabs = as.list(hzl),
  histlab = hl, citems = citems,
  bases = I(names(bm)), opacity = opac, def_id = dflt$id,
  extent = "CRA extent", about = about, repo = repo_url, cat = cat_url,
  base_credit = paste("Basemaps: Esri World Topographic Map and World Imagery",
                      "(© Esri and its data providers)."),
  bounds = c(list(list(g = "CRA extent", col = "#555555", dash = "4 4", on = FALSE)),
             lapply(bnd, function(b) b[c("g", "col", "dash", "on")])),
  aea = aea,
  data = as.character(toJSON(list(nx = nx, qd = qd, bnds = bnds, cmp = cmp),
                             auto_unbox = TRUE, na = "null", digits = NA, null = "null"))
)

css <- r"---[
.cra-panel{background:#fff;padding:0;border-radius:6px;box-shadow:0 1px 6px rgba(0,0,0,.35);
  font:12px/1.35 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;width:320px;max-width:80vw;color:#222}
.cra-hd{display:flex;justify-content:space-between;align-items:center;padding:7px 10px;
  font-weight:600;font-size:13px;border-bottom:1px solid #e3e3e3}
.cra-tg{border:0;background:none;font-size:16px;cursor:pointer;line-height:1;padding:0 4px}
.cra-bd{max-height:calc(100vh - 120px);overflow-y:auto;padding:6px 10px 8px}
.cra-sec{margin:4px 0 8px}
.cra-h{font-weight:600;text-transform:uppercase;letter-spacing:.04em;font-size:11px;color:#666;margin:6px 0 4px}
.cra-row{margin:3px 0;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
.cra-row select{max-width:210px;font-size:12px}
.cra-row label.dis{color:#aaa}
.cra-row input[type=range]{width:140px;vertical-align:middle}
.cra-sr{font-size:11px;color:#555;margin:6px 0 2px;font-style:italic}
.cra-panel details{margin:2px 0}
.cra-panel summary{cursor:pointer;padding:2px 0}
.cra-c>summary{font-weight:600}
.cra-c>details{margin-left:12px}
.cra-f>summary{color:#333}
.cra-ind{display:block;margin:2px 0 2px 16px;cursor:pointer}
.cra-ind.off{color:#aaa;cursor:not-allowed}
.cra-tag{font-size:10px;color:#777;border:1px solid #ccc;border-radius:3px;padding:0 3px;margin-left:3px}
.cra-sw{display:inline-block;width:22px;height:0;border-top-width:3px;border-top-style:solid;
  vertical-align:middle;margin-right:5px}
.cra-legend{background:#fff;padding:7px 10px;border-radius:6px;box-shadow:0 1px 6px rgba(0,0,0,.35);
  font:12px/1.3 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;width:260px}
.cra-legend .t{font-weight:600}
.cra-legend .s{color:#555;margin:2px 0 5px}
.cra-legend .bar{height:12px;border-radius:2px}
.cra-legend .ax{display:flex;justify-content:space-between;margin-top:2px}
.cra-legend .n{color:#777;font-size:11px;margin-top:3px}
.cra-legend .src{color:#555;font-size:11px;margin-top:5px;border-top:1px solid #eee;padding-top:4px}
.cra-legend .hk{display:flex;align-items:center;gap:6px;margin:3px 0;font-size:11px}
.cra-legend .hk canvas{border:1px solid #999;flex:none}
.cra-ab p{margin:3px 0 6px}
.cra-ab ul{margin:2px 0 6px;padding-left:16px}
.cra-ab li{margin:0 0 5px}
.cra-ab .m{color:#666;font-size:11px}
.cra-panel a,.cra-legend a,.cra-q a{color:#1f6fb4}
.cra-qbtn{background:#fff;width:34px;height:34px;border-radius:4px;border:2px solid rgba(0,0,0,.2);
  cursor:pointer;display:flex;align-items:center;justify-content:center;padding:0}
.cra-qbtn.on{background:#1f6fb4;border-color:#1f6fb4}
.cra-qbtn svg{width:20px;height:20px}
.leaflet-container.cra-qmode{cursor:crosshair}
.cra-q{background:#fff;border-radius:6px;box-shadow:0 1px 6px rgba(0,0,0,.35);
  font:12px/1.35 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;width:360px;max-width:85vw;color:#222}
.cra-q .qb{max-height:max(220px, calc(100vh - 470px));overflow-y:auto;padding:6px 10px 8px}
.cra-q table{border-collapse:collapse;width:100%;margin:2px 0 6px}
.cra-q td,.cra-q th{padding:2px 4px;border-bottom:1px solid #eee;text-align:left;vertical-align:top}
.cra-q th{font-weight:600;color:#555;font-size:11px}
.cra-q td.n{text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}
.cra-q .rg{color:#777;font-size:11px}
.cra-q .cl1{color:#1a7f37}.cra-q .cl2{color:#777}.cra-q .cl3{color:#b35806}
.cra-q .m{color:#666;font-size:11px}
.cra-q select{font-size:12px;max-width:220px}
.cra-q .chart text{font-size:10px;fill:#444}
]---"

js <- r"---[
function(el, x, cfg) {
  var map = this, lm = map.layerManager, idx = {}, D = JSON.parse(cfg.data);
  var N = D.nx, Q = D.qd, B = D.bnds || {}, C = D.cmp;
  var A = function(v) { return v == null ? [] : [].concat(v); };
  var G = function(g) { return lm.getLayerGroup(g, true); };
  var on = function(g, s) {
    var l = G(g);
    if (s && !map.hasLayer(l)) map.addLayer(l);
    if (!s && map.hasLayer(l)) map.removeLayer(l);
  };
  var each = function(sel, f) { Array.prototype.forEach.call(pd.querySelectorAll(sel), f); };
  var esc = function(s) {
    return String(s == null ? '' : s).replace(/[&<>"]/g, function(c) {
      return {'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}[c];
    });
  };
  var fmt = function(v) {
    if (v == null || !isFinite(v)) return '&ndash;';
    var a = Math.abs(v);
    return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) :
      a === 0 ? '0' : String(Number(v.toPrecision(2)));
  };
  var sfmt = function(v) { return v == null || !isFinite(v) ? '&ndash;' : v.toFixed(2); };
  var ok = function(v) { return v != null && isFinite(v); };
  var q7 = function(a, p) {
    var b = a.filter(ok).sort(function(x, y) { return x - y; }), n = b.length;
    if (!n) return null;
    var h = (n - 1) * p, lo = Math.floor(h);
    return b[lo] + (h - lo) * (b[Math.min(lo + 1, n - 1)] - b[lo]);
  };
  var link = function(u, t) { return '<a href="' + esc(u) + '" target="_blank" rel="noopener">' + esc(t || u) + '</a>'; };

  // ---- projection, grids and scoring ----
  var AE = (function() {
    var a = 6378137, f = 1 / 298.257223563, e2 = 2 * f - f * f, e = Math.sqrt(e2), r = Math.PI / 180;
    var q = function(s) { return (1 - e2) * (s / (1 - e2 * s * s) - (1 / (2 * e)) * Math.log((1 - e * s) / (1 + e * s))); };
    var mm = function(p) { var s = Math.sin(p); return Math.cos(p) / Math.sqrt(1 - e2 * s * s); };
    var p1 = cfg.aea.lat_1 * r, p2 = cfg.aea.lat_2 * r, p0 = cfg.aea.lat_0 * r, l0 = cfg.aea.lon_0 * r;
    var m1 = mm(p1), m2 = mm(p2), q1 = q(Math.sin(p1)), q2 = q(Math.sin(p2)), q0 = q(Math.sin(p0));
    var n = (m1 * m1 - m2 * m2) / (q2 - q1), Cc = m1 * m1 + n * q1, r0 = a * Math.sqrt(Cc - n * q0) / n;
    return {
      fwd: function(lon, lat) {
        var rho = a * Math.sqrt(Cc - n * q(Math.sin(lat * r))) / n, th = n * (lon * r - l0);
        return [rho * Math.sin(th), r0 - rho * Math.cos(th)];
      },
      inv: function(x, y) {
        var sg = n < 0 ? -1 : 1, rho = sg * Math.sqrt(x * x + (r0 - y) * (r0 - y));
        var th = Math.atan2(sg * x, sg * (r0 - y)), qq = (Cc - rho * rho * n * n / (a * a)) / n;
        var p = Math.asin(qq / 2);
        for (var k = 0; k < 20; k++) {
          var s = Math.sin(p), c = Math.cos(p), w = 1 - e2 * s * s;
          var dp = w * w / (2 * c) * (qq / (1 - e2) - s / w + (1 / (2 * e)) * Math.log((1 - e * s) / (1 + e * s)));
          p += dp;
          if (Math.abs(dp) < 1e-12) break;
        }
        return [(l0 + th / n) / r, p / r];
      }
    };
  })();
  var dec = {};
  var arr = function(id) {
    if (dec[id]) return dec[id];
    var z = Q.st[id], bin = atob(z.v), by = new Uint8Array(bin.length), s = null;
    for (var i = 0; i < by.length; i++) by[i] = bin.charCodeAt(i);
    if (z.s) { var b2 = atob(z.s); s = new Uint8Array(b2.length); for (var j = 0; j < s.length; j++) s[j] = b2.charCodeAt(j); }
    return (dec[id] = {u: new Float32Array(by.buffer), s: s});
  };
  var cell = function(ll) {
    var xy = AE.fwd(ll.lng, ll.lat), c = Math.floor((xy[0] - Q.x0) / Q.res), r = Math.floor((Q.y1 - xy[1]) / Q.res);
    if (c < 0 || r < 0 || c >= Q.nc || r >= Q.nr) return null;
    var ctr = AE.inv(Q.x0 + (c + 0.5) * Q.res, Q.y1 - (r + 0.5) * Q.res);
    return {k: r * Q.nc + c, lon: ctr[0], lat: ctr[1]};
  };
  var stv = function(id, k) {
    var u = arr(id).u[k];
    return isNaN(u) ? null : u;
  };
  var score = function(id, v) {
    var b = B[id];
    if (!ok(v) || !b || b.normalise !== 'minmax') return null;
    var t = b.transform === 'log1p' ? Math.log1p(v) : b.transform === 'sqrt' ? Math.sqrt(v) : v;
    var s = Math.min(1, Math.max(0, (t - b.lo) / (b.hi - b.lo)));
    return b.direction == -1 ? 1 - s : s;
  };
  var sts = function(id, k) {
    var a = arr(id);
    if (a.s) return a.s[k] === 255 ? null : a.s[k] / 250;
    return score(id, stv(id, k));
  };
  var bil = function(v, off, lon, lat) {
    var px = (lon - N.x0) / N.dx - 0.5, py = (N.y1 - lat) / N.dy - 0.5;
    var ic = Math.floor(px + 0.5), jc = Math.floor(py + 0.5);
    if (ic < 0 || jc < 0 || ic >= N.nc || jc >= N.nr || !ok(v[off + jc * N.nc + ic])) return null;
    var i0 = Math.floor(px), j0 = Math.floor(py), fx = px - i0, fy = py - j0, s = 0, w = 0;
    [[0, 0, (1 - fx) * (1 - fy)], [1, 0, fx * (1 - fy)], [0, 1, (1 - fx) * fy], [1, 1, fx * fy]].forEach(function(t) {
      var i = i0 + t[0], j = j0 + t[1];
      if (t[2] <= 0 || i < 0 || j < 0 || i >= N.nc || j >= N.nr) return;
      var z = v[off + j * N.nc + i];
      if (!ok(z)) return;
      s += z * t[2]; w += t[2];
    });
    return w > 0 ? s / w : null;
  };
  var nm = N ? N.models.length : 0, ncell = N ? N.nc * N.nr : 0;
  var gcmv = function(id, p, lon, lat) {
    var v = N.ind[id].v[p], o = [];
    for (var m = 0; m < nm; m++) o.push(bil(v, m * ncell, lon, lat));
    return o;
  };
  var mcache = {};
  var nmed = function(id, p) {
    var k = id + '|' + p;
    if (mcache[k]) return mcache[k];
    var v = N.ind[id].v[p], o = [];
    for (var c = 0; c < ncell; c++) {
      var z = [];
      for (var m = 0; m < nm; m++) z.push(v[m * ncell + c]);
      o.push(q7(z, 0.5));
    }
    return (mcache[k] = o);
  };
  var wmean = function(vals, w) {
    var num = 0, den = 0, tot = 0;
    for (var i = 0; i < vals.length; i++) { tot += w[i]; if (ok(vals[i])) { num += w[i] * vals[i]; den += w[i]; } }
    return den >= C.min_w * tot && den > 0 ? num / den : null;
  };
  var rfun = function(h, v, e) {
    if (!ok(h) || !ok(v) || !ok(e)) return null;
    return C.method === 'product' ? h * v * e : C.method === 'mean' ? (h + v + e) / 3 : Math.cbrt(h * v * e);
  };

  // ---- state and layer index ----
  var st = {id: cfg.def_id, kind: 'ind', view: 'value', ssp: A(cfg.ssps)[0], hz: 'hist', stat: 'median',
            op: cfg.opacity, hatch: true, cur: null, q: null, chartv: null};
  var hid = function() { return decodeURIComponent(location.hash.replace(/^#/, '')); };
  A(cfg.tree).forEach(function(sr) { A(sr.comps).forEach(function(c) { A(c.factors).forEach(function(f) {
    A(f.inds).forEach(function(i) { if (i.built) idx[i.id] = i; if (i.built && i.id === hid()) st.id = i.id; });
  }); }); });
  var cidx = {};
  A(cfg.citems).forEach(function(c) { cidx['cmp:' + c.key] = c; });
  if (cidx[hid()]) { st.id = hid(); st.kind = 'cmp'; }
  var period = function() {
    if (st.hz === 'hist') return cfg.pk.hist;
    return cfg.pk[st.hz + '|' + st.ssp];
  };

  A(cfg.layers).forEach(function(g) { on(g, false); });

  // ---- panel ----
  var opt = function(v, t, sel) { return '<option value="' + esc(v) + '"' + (sel ? ' selected' : '') + '>' + esc(t) + '</option>'; };
  var h = '<div class="cra-hd"><span>' + esc(cfg.title) +
    '</span><button type="button" class="cra-tg" title="Collapse">&#8722;</button></div><div class="cra-bd">';
  h += '<div class="cra-sec"><div class="cra-h">Display</div>' +
    '<div class="cra-row cra-vw"><label><input type="radio" name="cra-ty" value="value" checked> Values</label>' +
    '<label><input type="radio" name="cra-ty" value="score"> Scores (0&ndash;1)</label>' +
    '<label title="Median change of the GCMs from their own historical run, with IPCC AR6 robustness hatching">' +
    '<input type="radio" name="cra-ty" value="change"> Change from historical</label></div>' +
    '<div class="cra-row cra-strow" style="display:none"><label>Ensemble <select class="cra-stat"></select></label></div>' +
    '<div class="cra-row"><label>Period <select class="cra-hz">' + opt('hist', cfg.histlab, true) +
    A(cfg.hzs).map(function(z) { return opt(z, cfg.hzlabs[z], false); }).join('') + '</select></label></div>' +
    '<div class="cra-row"><label>Scenario <select class="cra-ssp">' +
    A(cfg.ssps).map(function(s, j) { return opt(s, cfg.ssplabs[s], j === 0); }).join('') + '</select></label></div>' +
    '<div class="cra-row cra-htrow"><label><input type="checkbox" class="cra-hatch" checked> Model agreement hatching</label></div>' +
    '<div class="cra-row"><label>Opacity <input type="range" class="cra-op" min="0" max="1" step="0.05" value="' +
    st.op + '"></label></div>' +
    '<div class="cra-row m" style="color:#666">Click the map (query tool on, top left) for values at a location.</div></div>';
  h += '<div class="cra-sec"><div class="cra-h">Layer</div>';
  h += '<label class="cra-ind" style="margin-left:0"><input type="radio" name="cra-ind" value=""> None</label>';
  if (A(cfg.citems).length) {
    h += '<details class="cra-c" open><summary>Composite indices (GCM ensemble)</summary>';
    A(cfg.citems).forEach(function(c) {
      var id = 'cmp:' + c.key;
      h += '<label class="cra-ind"><input type="radio" name="cra-ind" value="' + esc(id) + '"' +
        (id === st.id ? ' checked' : '') + '> ' + esc(c.label) + '</label>';
    });
    h += '</details>';
  }
  A(cfg.tree).forEach(function(sr) {
    h += '<div class="cra-sr" title="' + esc(sr.desc) + '">' + esc(sr.code) +
      (sr.desc ? ': ' + esc(sr.desc) : '') + '</div>';
    A(sr.comps).forEach(function(c) {
      h += '<details class="cra-c" open><summary>' + esc(c.name) + '</summary>';
      A(c.factors).forEach(function(f) {
        var fi = A(f.inds), op = fi.some(function(i) { return i.id === st.id; });
        h += '<details class="cra-f"' + (op ? ' open' : '') + '><summary>' + esc(f.name) + '</summary>';
        fi.forEach(function(i) {
          var tag = i.role && i.role !== 'Primary' ? '<span class="cra-tag">' + esc(String(i.role).toLowerCase()) + '</span>' : '';
          h += '<label class="cra-ind' + (i.built ? '' : ' off') + '" title="' +
            esc(i.built ? i.q : 'Not built yet (inactive in the indicator register)') + '">' +
            '<input type="radio" name="cra-ind" value="' + esc(i.id) + '"' + (i.built ? '' : ' disabled') +
            (i.id === st.id ? ' checked' : '') + '> ' + esc(i.label) + tag + '</label>';
        });
        h += '</details>';
      });
      h += '</details>';
    });
  });
  h += '</div><div class="cra-sec"><div class="cra-h">Boundaries</div>';
  A(cfg.bounds).forEach(function(b) {
    h += '<label class="cra-ind" style="margin-left:0"><input type="checkbox" class="cra-b" value="' + esc(b.g) + '"' +
      (b.on ? ' checked' : '') + '> <span class="cra-sw" style="border-top-color:' + b.col +
      (b.dash ? ';border-top-style:dashed' : '') + '"></span>' + esc(b.g) + '</label>';
  });
  h += '</div><div class="cra-sec"><div class="cra-h">Basemap</div>';
  A(cfg.bases).forEach(function(b, j) {
    h += '<label class="cra-ind" style="margin-left:0"><input type="radio" name="cra-bm" value="' + esc(b) + '"' +
      (j === 0 ? ' checked' : '') + '> ' + esc(b) + '</label>';
  });
  var seen = {}, srcs = [];
  A(cfg.tree).forEach(function(sr) { A(sr.comps).forEach(function(c) { A(c.factors).forEach(function(f) {
    A(f.inds).forEach(function(i) {
      var s = i.src || {}, k = (s.ds || '') + '|' + (s.pv || '');
      if (!i.built || !s.ds) return;
      if (!seen[k]) { seen[k] = {s: s, ids: []}; srcs.push(seen[k]); }
      seen[k].ids.push(i.label);
    });
  }); }); });
  h += '</div><details class="cra-sec cra-ab"><summary class="cra-h" style="display:list-item">About &amp; data sources</summary>' +
    '<p>' + esc(cfg.about) + '</p>' +
    (C ? '<p class="m">Composite: component index = weighted mean of indicator scores; risk = ' +
      esc(C.method === 'geomean' ? 'geometric mean of hazard, vulnerability and exposure' : C.method) +
      ', computed per GCM. &ldquo;High risk&rdquo; = at or above the ' + Math.round(C.q_high * 100) +
      'th percentile of historical median risk (' + sfmt(C.thr) + ').</p>' : '') +
    '<p class="m">Hatching (IPCC AR6 advanced method): a GCM&rsquo;s change is significant when it exceeds ' +
    '1.645 &times; its interannual standard deviation &times; &radic;(1/n<sub>hist</sub> + 1/n<sub>future</sub>) ' +
    '(rainfall CV: permutation test). No hatching = robust signal (&ge;66% of GCMs significant and &ge;80% agree on the sign); ' +
    'single lines = no change or no robust signal (&lt;66% significant); crossed lines = conflicting signals (&ge;66% significant, &lt;80% agree).</p>' +
    '<p><b>' + link(cfg.cat, 'Indicator catalogue') +
    '</b>: definitions, data sources, processing and status of every layer.</p><ul>';
  srcs.forEach(function(x) {
    var s = x.s, u = A(s.url);
    h += '<li><b>' + esc(s.pv || '') + '</b>' + (s.pv ? ': ' : '') + esc(s.ds) + (s.yr ? ' (' + esc(s.yr) + ')' : '') +
      '<div class="m">Used for: ' + esc(x.ids.join('; ')) + '</div>' +
      (s.lic ? '<div class="m">Access / licence: ' + esc(s.lic) + '</div>' : '') +
      (s.cite ? '<div class="m">Cite: ' + esc(s.cite) + '</div>' : '') +
      (u.length ? '<div class="m">' + u.map(function(v) { return link(v, 'Data link'); }).join(' &middot; ') + '</div>' : '') + '</li>';
  });
  h += '</ul><p class="m">' + esc(cfg.base_credit) + '</p><p class="m">Code and methods: ' + link(cfg.repo) + '</p></details></div>';

  var panel = L.control({position: 'topright'});
  panel.onAdd = function() {
    var d = L.DomUtil.create('div', 'cra-panel');
    d.innerHTML = h;
    L.DomEvent.disableClickPropagation(d);
    L.DomEvent.disableScrollPropagation(d);
    return d;
  };
  panel.addTo(map);
  var pd = panel.getContainer(), hzSel = pd.querySelector('.cra-hz'), sspSel = pd.querySelector('.cra-ssp'),
      statSel = pd.querySelector('.cra-stat');

  var legend = L.control({position: 'bottomleft'});
  legend.onAdd = function() { return L.DomUtil.create('div', 'cra-legend'); };
  legend.addTo(map);
  var ld = legend.getContainer();

  // ---- change layer (native grid, drawn client-side) ----
  map.createPane('craChg');
  map.getPane('craChg').style.zIndex = 390;
  var chg = {vals: null, cls: null, max: 1, worse: 1};
  var DIV = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#f7f7f7', '#fddbc7', '#f4a582', '#d6604d', '#b2182b'];
  var hex = function(c) { return [1, 3, 5].map(function(i) { return parseInt(c.substr(i, 2), 16); }); };
  var divcol = function(t) {
    t = Math.max(0, Math.min(1, t)) * (DIV.length - 1);
    var i = Math.min(Math.floor(t), DIV.length - 2), f = t - i, a = hex(DIV[i]), b = hex(DIV[i + 1]);
    return 'rgb(' + a.map(function(v, k) { return Math.round(v + f * (b[k] - v)); }).join(',') + ')';
  };
  var chgCol = function(v) { return divcol(0.5 + 0.5 * chg.worse * v / chg.max); };
  var hatch = function(ctx, x0, y0, w, hh, ox, oy, cl) {
    var sp = 7;
    ctx.save(); ctx.beginPath(); ctx.rect(x0, y0, w, hh); ctx.clip();
    ctx.strokeStyle = 'rgba(0,0,0,0.55)'; ctx.lineWidth = 1; ctx.beginPath();
    var o1 = (((-ox - oy) % sp) + sp) % sp;
    for (var t = o1; t <= 512; t += sp) { ctx.moveTo(t, 0); ctx.lineTo(t - 256, 256); }
    if (cl === 3) {
      var o2 = (((oy - ox) % sp) + sp) % sp;
      for (var u = o2 - 256; u <= 256; u += sp) { ctx.moveTo(u, 0); ctx.lineTo(u + 256, 256); }
    }
    ctx.stroke(); ctx.restore();
  };
  var ChgLayer = L.GridLayer.extend({
    createTile: function(c) {
      var t = L.DomUtil.create('canvas'), sz = this.getTileSize();
      t.width = sz.x; t.height = sz.y;
      if (!chg.vals) return t;
      var ctx = t.getContext('2d'), z = c.z, ox = c.x * sz.x, oy = c.y * sz.y;
      for (var j = 0; j < N.nr; j++) for (var i = 0; i < N.nc; i++) {
        var k = j * N.nc + i, v = chg.vals[k];
        if (!ok(v)) continue;
        var p1 = map.project([N.y1 - j * N.dy, N.x0 + i * N.dx], z), p2 = map.project([N.y1 - (j + 1) * N.dy, N.x0 + (i + 1) * N.dx], z);
        var x0 = p1.x - ox, y0 = p1.y - oy, w = p2.x - p1.x, hh = p2.y - p1.y;
        if (x0 > sz.x || y0 > sz.y || x0 + w < 0 || y0 + hh < 0) continue;
        ctx.globalAlpha = st.op; ctx.fillStyle = chgCol(v); ctx.fillRect(x0, y0, w, hh);
        ctx.globalAlpha = 1;
        var cl = chg.cls ? chg.cls[k] : 0;
        if (st.hatch && (cl === 2 || cl === 3)) hatch(ctx, x0, y0, w, hh, ox, oy, cl);
      }
      return t;
    }
  });
  var chgL = new ChgLayer({pane: 'craChg', tileSize: 256});
  var chgOn = false;
  var swatch = function(cl) {
    var c = document.createElement('canvas'); c.width = 26; c.height = 14;
    var x = c.getContext('2d'); x.fillStyle = '#f7f7f7'; x.fillRect(0, 0, 26, 14);
    if (cl > 1) hatch(x, 0, 0, 26, 14, 0, 0, cl);
    return c.toDataURL();
  };
  var SW = {1: swatch(1), 2: swatch(2), 3: swatch(3)};
  var CLS = {1: 'Robust signal', 2: 'No change or no robust signal', 3: 'Conflicting signals'};

  // per-GCM change at native cells for indicator id, period p
  var nchg = function(id, p) {
    var d = N.ind[id], fv = d.v[p], hv = d.v[cfg.pk.hist], out = [];
    for (var k = 0; k < ncell; k++) {
      var dd = [];
      for (var m = 0; m < nm; m++) {
        var a = fv[m * ncell + k], b = hv[m * ncell + k];
        if (ok(a) && ok(b)) dd.push(d.rel ? (b > 0 ? 100 * (a - b) / b : null) : a - b);
      }
      out.push(q7(dd, 0.5));
    }
    return out;
  };
  var chgMax = {};
  var maxAbs = function(id) {
    if (chgMax[id]) return chgMax[id];
    var mx = 0;
    A(cfg.periods).forEach(function(p) {
      if (p === cfg.pk.hist || !N.ind[id].v[p]) return;
      nchg(id, p).forEach(function(v) { if (ok(v)) mx = Math.max(mx, Math.abs(v)); });
    });
    return (chgMax[id] = mx || 1);
  };

  var front = function() { if (st.cur) G(st.cur).eachLayer(function(l) { if (l.bringToFront) l.bringToFront(); }); };
  var setOp = function() {
    if (st.cur) G(st.cur).eachLayer(function(l) { if (l.setOpacity) l.setOpacity(st.op); });
    if (chgOn) chgL.redraw();
  };
  var srcLine = function(i) {
    var s = i.src || {}, u = A(s.url);
    if (!s.ds) return '';
    var t = 'Source: ' + (s.pv ? s.pv + ', ' : '') + s.ds;
    return '<div class="src">' + (u.length ? link(u[0], t) : esc(t)) + (s.lic ? '<br>' + esc(s.lic) : '') + '</div>';
  };
  var catLine = function(i) {
    return '<div class="src"><a href="' + esc(cfg.cat) + '#' + encodeURIComponent(i.id) +
      '" target="_blank" rel="noopener">Layer details in the indicator catalogue &rarr;</a></div>';
  };
  if (map.attributionControl) map.attributionControl.addAttribution(
    'CRA layers: ' + link(cfg.repo, 'Umzimvubu CRA') + ' (sources under About)');

  var setDis = function(sel, d) { sel.disabled = d; sel.parentNode.className = d ? 'dis' : ''; };

  function draw() {
    if (st.cur) on(st.cur, false);
    st.cur = null;
    if (chgOn) { map.removeLayer(chgL); chgOn = false; }
    var vw = pd.querySelector('.cra-vw'), sr = pd.querySelector('.cra-strow'), hr = pd.querySelector('.cra-htrow');
    var i = st.kind === 'cmp' ? cidx[st.id] : idx[st.id];
    if (!i) { ld.style.display = 'none'; setDis(hzSel, true); setDis(sspSel, true); renderQ(); return; }
    var ps = A(i.periods), stat = st.kind === 'cmp';
    var dynamic = ps.indexOf('baseline') < 0;
    var canChg = !stat && i.nex && N && N.ind[i.id];
    vw.style.display = stat ? 'none' : '';
    sr.style.display = stat ? '' : 'none';
    each('input[name=cra-ty]', function(e) {
      if (e.value === 'change') { e.disabled = !canChg; e.parentNode.className = canChg ? '' : 'dis'; }
    });
    if (st.view === 'change' && !canChg) {
      st.view = 'value';
      each('input[name=cra-ty]', function(e) { e.checked = e.value === 'value'; });
    }
    var isChg = !stat && st.view === 'change';
    if (isChg && st.hz === 'hist') { st.hz = A(cfg.hzs)[A(cfg.hzs).length - 1]; hzSel.value = st.hz; }
    hzSel.options[0].disabled = isChg;
    setDis(hzSel, !dynamic);
    setDis(sspSel, !dynamic || st.hz === 'hist');
    hr.style.display = isChg ? '' : 'none';
    if (stat) {
      if (A(i.stats).indexOf(st.stat) < 0) st.stat = A(i.stats)[0];
      statSel.innerHTML = A(i.stats).map(function(s) { return opt(s, i.slabs[s], s === st.stat); }).join('');
    }
    var p = dynamic ? period() : 'baseline';
    if (ps.indexOf(p) < 0) p = ps[0];
    ld.style.display = '';
    if (isChg) {
      var d = N.ind[i.id];
      chg.vals = nchg(i.id, p); chg.cls = d.cls[p] ? A(d.cls[p]) : null;
      chg.max = maxAbs(i.id); chg.worse = i.dir == -1 ? -1 : 1;
      map.addLayer(chgL); chgOn = true; chgL.redraw();
      var u = d.rel ? '%' : (i.units || '');
      var cols = []; for (var c = 0; c <= 8; c++) cols.push(chgCol(chg.max * (c / 4 - 1)));
      var hk = chg.cls ? [1, 2, 3].map(function(k) {
        return '<div class="hk"><img src="' + SW[k] + '" width="26" height="14" style="border:1px solid #999"> ' + CLS[k] + '</div>';
      }).join('') : '<div class="n">Robustness hatching not available (re-run 201).</div>';
      ld.innerHTML = '<div class="t">' + esc(i.label) + '</div><div class="s">Median change of ' + nm + ' GCMs from historical' +
        (u ? ' (' + esc(u) + ')' : '') + '<br>' + esc(cfg.plabs[p] || p) + '</div>' +
        '<div class="bar" style="background:linear-gradient(to right,' + cols.join(',') + ')"></div>' +
        '<div class="ax"><span>' + fmt(-chg.max) + '</span><span>0</span><span>+' + fmt(chg.max) + '</span></div>' +
        '<div class="n">Red = change towards higher risk. Native ~25 km climate grid.</div>' + hk +
        '<div class="n">IPCC AR6 method; see About.</div>' + srcLine(i) + catLine(i);
    } else {
      var key = stat ? ['cmp', i.key, st.stat].join('|') : [st.id, st.view].join('|');
      st.cur = key + '|' + p;
      on(st.cur, true);
      setOp();
      front();
      var lg = cfg.legends[key];
      var sub = stat ? (st.stat === 'high_share' ? '% of GCMs at or above the high-risk threshold' :
                        st.stat === 'index' ? 'Index: 0 = lowest, 1 = highest' : i.slabs[st.stat] + ' (0 = lowest, 1 = highest)') :
        st.view === 'score' ? 'Score: 0 = lowest, 1 = highest risk' : 'Value' + (i.units ? ' (' + i.units + ')' : '') +
        (i.nex ? ': median of ' + nm + ' GCMs' : '');
      var note = !stat && st.view === 'value' ? (i.dir == -1 ? 'Lower values = higher risk' : 'Higher values = higher risk') : '';
      ld.innerHTML = '<div class="t">' + esc(stat ? i.label : i.label) + '</div><div class="s">' + esc(sub) + '<br>' +
        esc(cfg.plabs[p] || p) + '</div><div class="bar" style="background:linear-gradient(to right,' +
        A(lg.cols).join(',') + ')"></div><div class="ax"><span>' + fmt(lg.min) + '</span><span>' +
        fmt(lg.max) + '</span></div>' + (note ? '<div class="n">' + note + '</div>' : '') +
        (stat ? '' : srcLine(i) + catLine(i));
    }
    renderQ();
  }

  // ---- location query ----
  var qmode = true, qmk = null;
  var qbtn = L.control({position: 'topleft'});
  qbtn.onAdd = function() {
    var b = L.DomUtil.create('button', 'cra-qbtn on');
    b.type = 'button'; b.title = 'Query a location (click the map)';
    b.innerHTML = '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="7"/>' +
      '<path d="M12 2v5M12 17v5M2 12h5M17 12h5"/></svg>';
    L.DomEvent.disableClickPropagation(b);
    L.DomEvent.on(b, 'click', function() {
      qmode = !qmode;
      b.className = 'cra-qbtn' + (qmode ? ' on' : '');
      b.style.color = qmode ? '#fff' : '#333';
      map.getContainer().classList.toggle('cra-qmode', qmode);
    });
    b.style.color = '#fff';
    return b;
  };
  qbtn.addTo(map);
  map.getContainer().classList.add('cra-qmode');
  var qc = L.control({position: 'topleft'});
  qc.onAdd = function() {
    var d = L.DomUtil.create('div', 'cra-q');
    d.style.display = 'none';
    L.DomEvent.disableClickPropagation(d);
    L.DomEvent.disableScrollPropagation(d);
    return d;
  };
  qc.addTo(map);
  var qd = qc.getContainer();

  var compute = function(cl, p) {
    var o = {ind: {}, cmp: null};
    Object.keys(idx).forEach(function(id) {
      var i = idx[id];
      if (i.nex && N && N.ind[id] && N.ind[id].v[p]) {
        var g = gcmv(id, p, cl.lon, cl.lat), s = g.map(function(v) { return score(id, v); });
        var med = bil(nmed(id, p), 0, cl.lon, cl.lat);
        o.ind[id] = {g: g, s: s, v: med, p10: q7(g, 0.1), p90: q7(g, 0.9), sc: score(id, med)};
      } else if (!i.nex && Q.st[id]) {
        var v = stv(id, cl.k);
        o.ind[id] = {v: v, sc: sts(id, cl.k)};
      }
    });
    if (C) {
      var cv = {}, comps = C.comps, hz = comps.hazard;
      ['vulnerability', 'exposure'].forEach(function(cp) {
        var z = comps[cp];
        cv[cp] = z ? wmean(A(z.ids).map(function(id) { return o.ind[id] ? o.ind[id].sc : null; }), A(z.w)) : null;
      });
      var hm = [], rm = [];
      for (var m = 0; m < nm; m++) {
        var hv = wmean(A(hz.ids).map(function(id) {
          var x = o.ind[id];
          return !x ? null : x.s ? x.s[m] : x.sc;
        }), A(hz.w));
        hm.push(hv); rm.push(rfun(hv, cv.vulnerability, cv.exposure));
      }
      var rv = rm.filter(ok);
      o.cmp = {h: hm, r: rm, v: cv.vulnerability, e: cv.exposure,
               hq: [q7(hm, 0.5), q7(hm, 0.1), q7(hm, 0.9)], rq: [q7(rm, 0.5), q7(rm, 0.1), q7(rm, 0.9)],
               hs: rv.length ? 100 * rv.filter(function(v) { return v >= C.thr; }).length / rv.length : null};
    }
    return o;
  };
  var gcmChg = function(id, p, cl) {
    var d = N.ind[id], f = gcmv(id, p, cl.lon, cl.lat), b = gcmv(id, cfg.pk.hist, cl.lon, cl.lat);
    return f.map(function(v, m) { return ok(v) && ok(b[m]) ? (d.rel ? (b[m] > 0 ? 100 * (v - b[m]) / b[m] : null) : v - b[m]) : null; });
  };
  var ncellAt = function(lon, lat) {
    var i = Math.floor((lon - N.x0) / N.dx), j = Math.floor((N.y1 - lat) / N.dy);
    return i < 0 || j < 0 || i >= N.nc || j >= N.nr ? -1 : j * N.nc + i;
  };
  var chart = function(rows, xl, thr) {
    var W = 330, rh = 22, top = 6, lw = 70, H = top + rows.length * rh + 26, all = [];
    rows.forEach(function(r) { r.v.forEach(function(v) { if (ok(v)) all.push(v); }); });
    if (thr != null) all.push(thr);
    if (!all.length) return '<div class="m">No GCM values at this location.</div>';
    var lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
    if (hi === lo) { lo -= 0.5; hi += 0.5; }
    var pad = (hi - lo) * 0.05; lo -= pad; hi += pad;
    var X = function(v) { return lw + (v - lo) / (hi - lo) * (W - lw - 8); };
    var s = '<svg class="chart" width="' + W + '" height="' + H + '" viewBox="0 0 ' + W + ' ' + H + '">';
    if (thr != null) s += '<line x1="' + X(thr) + '" x2="' + X(thr) + '" y1="' + top + '" y2="' + (top + rows.length * rh) +
      '" stroke="#b2182b" stroke-dasharray="3 3"><title>High-risk threshold ' + sfmt(thr) + '</title></line>';
    if (lo < 0 && hi > 0 && thr == null) s += '<line x1="' + X(0) + '" x2="' + X(0) + '" y1="' + top + '" y2="' + (top + rows.length * rh) + '" stroke="#999"/>';
    rows.forEach(function(r, j) {
      var y = top + j * rh + rh / 2, v = r.v.filter(ok);
      s += '<text x="0" y="' + (y + 3) + '">' + esc(r.lab) + '</text>';
      if (!v.length) return;
      var a = q7(v, 0.1), b = q7(v, 0.9), md = q7(v, 0.5);
      s += '<line x1="' + X(a) + '" x2="' + X(b) + '" y1="' + y + '" y2="' + y + '" stroke="#bbb" stroke-width="6" stroke-linecap="round"><title>P10&ndash;P90: ' +
        fmt(a) + ' to ' + fmt(b) + '</title></line>';
      r.v.forEach(function(z, m) {
        if (!ok(z)) return;
        var f = r.f ? r.f[m] : true;
        s += '<circle cx="' + X(z).toFixed(1) + '" cy="' + (y + ((m % 3) - 1) * 3) + '" r="2.6" fill="' + (f ? r.col : '#fff') +
          '" stroke="' + r.col + '" stroke-width="1"><title>' + esc(N.models[m]) + ': ' + fmt(z) +
          (r.f ? (f ? ' (significant change)' : ' (not significant)') : '') + '</title></circle>';
      });
      s += '<line x1="' + X(md) + '" x2="' + X(md) + '" y1="' + (y - 7) + '" y2="' + (y + 7) + '" stroke="#111" stroke-width="2"><title>Median ' + fmt(md) + '</title></line>';
    });
    var yb = top + rows.length * rh + 4;
    s += '<line x1="' + lw + '" x2="' + (W - 8) + '" y1="' + yb + '" y2="' + yb + '" stroke="#999"/>';
    [lo + pad, (lo + hi) / 2, hi - pad].forEach(function(v) {
      s += '<text x="' + X(v) + '" y="' + (yb + 12) + '" text-anchor="middle">' + fmt(v) + '</text>';
    });
    s += '<text x="' + ((lw + W) / 2) + '" y="' + (yb + 22) + '" text-anchor="middle">' + esc(xl) + '</text></svg>';
    return s;
  };
  var SCOL = {hist: '#555', ssp126: '#1a9850', ssp245: '#2c7bb6', ssp370: '#f46d43', ssp585: '#a50026'};

  function renderQ() {
    if (!st.q) { qd.style.display = 'none'; return; }
    var cl = st.q, p = period() || cfg.pk.hist, hz = st.hz === 'hist' ? A(cfg.hzs)[A(cfg.hzs).length - 1] : st.hz;
    var o = compute(cl, p), s = '';
    s += '<div class="cra-hd"><span>Location ' + cl.lat.toFixed(4) + ', ' + cl.lon.toFixed(4) +
      '</span><button type="button" class="cra-tg cra-qx" title="Close">&times;</button></div><div class="qb">';
    s += '<div class="m">1 km cell centre; climate values interpolated from the ~25 km GCM grid as in the pipeline. ' +
      esc(cfg.plabs[p] || p) + '.</div>';
    if (o.cmp) {
      var c = o.cmp, rg = function(q) { return '<span class="rg">[' + sfmt(q[1]) + '&ndash;' + sfmt(q[2]) + ']</span>'; };
      s += '<div class="cra-h">Composite (median of ' + nm + ' GCMs [P10&ndash;P90])</div><table>' +
        '<tr><td><b>Climate risk</b></td><td class="n"><b>' + sfmt(c.rq[0]) + '</b> ' + rg(c.rq) + '</td></tr>' +
        '<tr><td>GCMs with high risk</td><td class="n">' + (ok(c.hs) ? Math.round(c.hs) + '%' : '&ndash;') + '</td></tr>' +
        '<tr><td>Hazard index</td><td class="n">' + sfmt(c.hq[0]) + ' ' + rg(c.hq) + '</td></tr>' +
        '<tr><td>Vulnerability index</td><td class="n">' + sfmt(c.v) + '</td></tr>' +
        '<tr><td>Exposure index</td><td class="n">' + sfmt(c.e) + '</td></tr></table>';
    }
    var nk = N ? ncellAt(cl.lon, cl.lat) : -1;
    if (N) {
      var opts = [];
      if (C) { opts.push(['risk', 'Climate risk']); opts.push(['hazard', 'Hazard index']); }
      Object.keys(idx).forEach(function(id) { if (idx[id].nex && N.ind[id]) { opts.push([id, idx[id].label]); opts.push([id + '|chg', idx[id].label + ' (change)']); } });
      if (!st.chartv || !opts.some(function(z) { return z[0] === st.chartv; }))
        st.chartv = st.kind === 'ind' && idx[st.id] && idx[st.id].nex ? st.id + (st.view === 'change' ? '|chg' : '') : opts[0][0];
      s += '<div class="cra-h">Spread across GCMs (' + esc(cfg.hzlabs[hz] || hz) + ')</div>' +
        '<div class="cra-row"><select class="cra-chv">' + opts.map(function(z) { return opt(z[0], z[1], z[0] === st.chartv); }).join('') + '</select></div>';
      var cv = st.chartv, isC = /\|chg$/.test(cv), id = cv.replace(/\|chg$/, ''), rows2 = [], thr = null, xl = '';
      var pl = [];
      if (!isC) pl.push(['hist', cfg.pk.hist]);
      A(cfg.ssps).forEach(function(sp) { var q = cfg.pk[hz + '|' + sp]; if (q) pl.push([sp, q]); });
      pl.forEach(function(z) {
        var lab = z[0] === 'hist' ? 'Historical' : cfg.ssplabs[z[0]], v, f = null;
        if (cv === 'risk' || cv === 'hazard') {
          var oo = compute(cl, z[1]);
          v = cv === 'risk' ? oo.cmp.r : oo.cmp.h;
          if (cv === 'risk') f = v.map(function(x) { return ok(x) && x >= C.thr; });
        } else if (isC) {
          v = gcmChg(id, z[1], cl);
          var sg = N.ind[id].sig[z[1]];
          if (sg && nk >= 0) f = N.models.map(function(m, j) { return sg.charAt(j * ncell + nk) === '1'; });
        } else {
          v = gcmv(id, z[1], cl.lon, cl.lat);
        }
        rows2.push({lab: lab, v: v, f: f, col: SCOL[z[0]] || '#333'});
      });
      var nt = '';
      if (cv === 'risk') { thr = C.thr; xl = 'Risk (0-1)'; nt = ' Dashed line = high-risk threshold; filled dot = GCM with high risk.'; }
      else if (cv === 'hazard') xl = 'Hazard index (0-1)';
      else if (isC) { xl = 'Change from historical' + (N.ind[id].rel ? ' (%)' : idx[id].units ? ' (' + idx[id].units + ')' : '');
        nt = ' Filled dot = change significant against that GCM&rsquo;s own variability (climate cell containing the point).'; }
      else xl = idx[id].units || idx[id].label;
      s += chart(rows2, xl, thr) + '<div class="m">Dots = GCMs, grey bar = P10&ndash;P90, black tick = median; hover for model names.' + nt + '</div>';
    }
    A(cfg.tree).forEach(function(sr) { A(sr.comps).forEach(function(cc) {
      var rows = '';
      A(cc.factors).forEach(function(f) { A(f.inds).forEach(function(i) {
        var x = o.ind[i.id];
        if (!i.built || !x) return;
        var u = i.units ? ' ' + esc(i.units) : '', extra = '';
        if (i.nex && st.hz !== 'hist') {
          var d = N.ind[i.id], ch = q7(gcmChg(i.id, p, cl), 0.5), k = d.cls[p] && nk >= 0 ? A(d.cls[p])[nk] : 0;
          extra = '<div class="rg">change ' + (ok(ch) && ch > 0 ? '+' : '') + fmt(ch) + (d.rel ? '%' : u) +
            (k ? ' &middot; <span class="cl' + k + '">' + CLS[k] + '</span>' : '') + '</div>';
        }
        rows += '<tr><td>' + esc(i.label) + extra + '</td><td class="n">' + fmt(x.v) + u +
          (i.nex ? '<div class="rg">[' + fmt(x.p10) + '&ndash;' + fmt(x.p90) + ']</div>' : '') +
          '</td><td class="n">' + sfmt(x.sc) + '</td></tr>';
      }); });
      if (rows) s += '<div class="cra-h">' + esc(cc.name) + '</div><table><tr><th>Indicator</th><th>Value</th><th>Score</th></tr>' + rows + '</table>';
    }); });
    s += '</div>';
    qd.innerHTML = s;
    qd.style.display = '';
    qd.querySelector('.cra-qx').addEventListener('click', function() {
      st.q = null; qd.style.display = 'none';
      if (qmk) { map.removeLayer(qmk); qmk = null; }
    });
    var cs = qd.querySelector('.cra-chv');
    if (cs) cs.addEventListener('change', function() { st.chartv = this.value; renderQ(); });
  }
  map.on('click', function(ev) {
    if (!qmode) return;
    var cl = cell(ev.latlng);
    if (!cl) return;
    st.q = cl;
    if (qmk) map.removeLayer(qmk);
    qmk = L.circleMarker([cl.lat, cl.lon], {radius: 6, color: '#111', weight: 2, fillColor: '#ffeb3b', fillOpacity: 1}).addTo(map);
    renderQ();
  });

  // ---- events ----
  each('input[name=cra-ind]', function(e) { e.addEventListener('change', function() {
    st.id = this.value; st.kind = /^cmp:/.test(st.id) ? 'cmp' : 'ind'; draw();
  }); });
  each('input[name=cra-ty]', function(e) { e.addEventListener('change', function() { st.view = this.value; draw(); }); });
  hzSel.addEventListener('change', function() { st.hz = this.value; draw(); });
  sspSel.addEventListener('change', function() { st.ssp = this.value; draw(); });
  statSel.addEventListener('change', function() { st.stat = this.value; draw(); });
  pd.querySelector('.cra-hatch').addEventListener('change', function() { st.hatch = this.checked; if (chgOn) chgL.redraw(); });
  pd.querySelector('.cra-op').addEventListener('input', function() { st.op = +this.value; setOp(); });
  each('input.cra-b', function(e) {
    on(e.value, e.checked);
    e.addEventListener('change', function() { on(this.value, this.checked); });
  });
  var base = function(b) {
    A(cfg.bases).forEach(function(x) { on(x, x === b); });
    G(b).eachLayer(function(l) { if (l.bringToBack) l.bringToBack(); });
    front();
  };
  each('input[name=cra-bm]', function(e) { e.addEventListener('change', function() { base(this.value); }); });
  var tg = pd.querySelector('.cra-tg'), bd = pd.querySelector('.cra-bd');
  tg.addEventListener('click', function() {
    var c = bd.style.display !== 'none';
    bd.style.display = c ? 'none' : '';
    tg.innerHTML = c ? '&#43;' : '&#8722;';
    tg.title = c ? 'Expand' : 'Collapse';
  });

  window.addEventListener('hashchange', function() {
    var k = hid();
    if (!idx[k] && !cidx[k]) return;
    st.id = k; st.kind = cidx[k] ? 'cmp' : 'ind';
    each('input[name=cra-ind]', function(e) {
      e.checked = e.value === k;
      if (e.checked && e.closest('details.cra-f')) e.closest('details.cra-f').open = true;
    });
    draw();
  });

  window.CRA = {st: st, compute: compute, cell: cell, AE: AE, draw: draw, period: period};
  base(A(cfg.bases)[0]);
  draw();
}
]---"

m <- m |>
  hideGroup(c(lyr, names(bm)[-1], "CRA extent")) |>
  fitBounds(roi[["xmin"]], roi[["ymin"]], roi[["xmax"]], roi[["ymax"]]) |>
  addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE)) |>
  prependContent(tags$style(HTML(css))) |>
  onRender(js, data = cfg)

saveWidget(m, out_html, selfcontained = FALSE, title = cfg$title)
