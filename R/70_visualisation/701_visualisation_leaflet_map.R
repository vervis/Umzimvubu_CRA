# ============================================================================
# 701_visualisation_leaflet_map.R
#
# Interactive Leaflet map of the stage-30 indicator layers
# (data/processed/indicators/<component>/<id>_value.tif and _score.tif),
# organised by the CRA framework hierarchy in CRA_data_sources_Claude.xlsx
# (sheet Indicators): Sub-risk > Component > Factor > Indicator.
#   1. Settings, boundaries and CRA extent.
#   2. Framework: xlsx hierarchy joined to the register on Register ID; rows
#      without stage-30 outputs are listed but disabled. Falls back to the
#      register (component/label) if the xlsx is missing.
#   3. Rasters: one hidden group per indicator x {value, score} x period.
#      Value palettes (viridis) are pooled across periods so change between
#      periods stays visible; scores use a fixed 0-1 YlOrRd scale.
#   4. Map: CARTO / Esri basemaps (no Referer needed, unlike the OSM tile
#      servers, which block file:// pages), catchment, district and
#      Matatiele boundaries as separate coloured groups, and a custom panel
#      (onRender JS): collapsible indicator tree with one active layer,
#      value/score switch, period selector, opacity slider, boundary and
#      basemap toggles, and a legend for the active layer.
#      Saved to outputs/maps/umzimvubu_cra_indicators_leaflet.html.
# ============================================================================

library(sf)
library(terra)
library(leaflet)
library(htmlwidgets)
library(htmltools)
library(readxl)

source("R/functions/pipeline_helpers.R")

reg_csv  <- "R/30_indicators/indicator_register.csv"
fw_xlsx  <- "../CRA_data_sources_Claude.xlsx"
ind_dir  <- "data/processed/indicators"
out_html <- "outputs/maps/umzimvubu_cra_indicators_leaflet.html"
opac     <- 0.8
dir.create(dirname(out_html), recursive = TRUE, showWarnings = FALSE)

rd <- function(f) st_transform(st_zm(st_make_valid(st_read(f, quiet = TRUE))), 4326)

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

if (file.exists(fw_xlsx)) {
  x  <- as.data.frame(read_excel(fw_xlsx, sheet = "Indicators", skip = 1))
  x  <- x[!is.na(x[["Register ID"]]), ]
  fw <- data.frame(id = x[["Register ID"]], subrisk = x[["Sub-risk"]],
                   component = x[["Component"]], factor = x[["Factor"]],
                   q = x[["Question answered"]], role = x[["Role"]])
  lk  <- as.data.frame(read_excel(fw_xlsx, sheet = "Lookups", col_names = FALSE,
                                  .name_repair = "minimal"))
  lk  <- lk[grepl("^SR\\d", lk[[1]]), ]
  srl <- setNames(as.character(lk[[2]]), lk[[1]])
} else {
  fw <- data.frame(id = reg$id, subrisk = reg$subrisk,
                   component = tools::toTitleCase(reg$component),
                   factor = reg$label, q = NA, role = NA)
  srl <- character()
}

k  <- match(fw$id, reg$id)
fw <- fw[!is.na(k), ]
k  <- k[!is.na(k)]
ind <- cbind(fw, reg[k, c("label", "units", "temporal", "direction")],
             cdir = tolower(reg$component[k]))
ind$label <- ifelse(is.na(ind$label), ind$id, ind$label)
ind$fv <- file.path(ind_dir, ind$cdir, paste0(ind$id, "_value.tif"))
ind$fs <- file.path(ind_dir, ind$cdir, paste0(ind$id, "_score.tif"))
ind$built <- file.exists(ind$fv) & file.exists(ind$fs)
if (!any(ind$built)) stop("No indicator outputs found in ", ind_dir)

plab <- function(p) {
  if (p == "baseline") return("Baseline (static layer)")
  y <- sub("_", "-", regmatches(p, regexpr("\\d{4}_\\d{4}", p)))
  s <- regmatches(p, regexpr("ssp\\d{3}", p))
  if (!length(s)) return(paste("Historical", y))
  h <- c(near = "Near", mid = "Mid", late = "Late")[sub("_.*", "", p)]
  sprintf("%s %s, SSP%s-%s.%s", h, y, substr(s, 4, 4), substr(s, 5, 5), substr(s, 6, 6))
}


# 3. Rasters + 4. Map ------------------------------------------------------------

bm <- c("Light (CARTO)" = "CartoDB.Positron",
        "Topographic (Esri)" = "Esri.WorldTopoMap",
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
for (i in which(ind$built)) {
  for (ty in c("value", "score")) {
    r  <- rast(ind[[if (ty == "value") "fv" else "fs"]][i])
    v  <- values(r, mat = FALSE)
    rg <- if (ty == "score") c(0, 1) else range(v[is.finite(v)])
    if (!all(is.finite(rg))) next
    if (diff(rg) == 0) rg <- rg + c(-0.5, 0.5)
    pal <- colorNumeric(if (ty == "score") "YlOrRd" else "viridis", rg,
                        na.color = "transparent")
    for (j in seq_len(nlyr(r))) {
      g <- paste(ind$id[i], ty, names(r)[j], sep = "|")
      m <- addRasterImage(m, r[[j]], colors = pal, opacity = opac,
                          project = TRUE, method = "ngb", group = g)
      lyr <- c(lyr, g)
    }
    leg[[paste(ind$id[i], ty, sep = "|")]] <-
      list(cols = pal(seq(rg[1], rg[2], length.out = 9)), min = rg[1], max = rg[2])
    per[[ind$id[i]]] <- names(r)
  }
}
ind$built <- ind$id %in% names(per)

pers <- unique(unlist(per))
pers <- c(setdiff(pers, "baseline"), intersect(pers, "baseline"))
dflt <- ind[ind$built, ][1, ]

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
             periods = I(if (d$built[j]) per[[d$id[j]]] else character())
           )))
         }))
       }))
})

cfg <- list(
  title = "Umzimvubu CRA indicators", tree = tree, layers = I(lyr), legends = leg,
  periods = I(pers), plabs = setNames(lapply(pers, plab), pers),
  bases = I(names(bm)), opacity = opac, def_id = dflt$id,
  def_period = per[[dflt$id]][1], extent = "CRA extent",
  bounds = c(list(list(g = "CRA extent", col = "#555555", dash = "4 4", on = FALSE)),
             lapply(bnd, function(b) b[c("g", "col", "dash", "on")]))
)

css <- r"---[
.cra-panel{background:#fff;padding:0;border-radius:6px;box-shadow:0 1px 6px rgba(0,0,0,.35);
  font:12px/1.35 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;width:310px;max-width:80vw;color:#222}
.cra-hd{display:flex;justify-content:space-between;align-items:center;padding:7px 10px;
  font-weight:600;font-size:13px;border-bottom:1px solid #e3e3e3}
.cra-tg{border:0;background:none;font-size:16px;cursor:pointer;line-height:1;padding:0 4px}
.cra-bd{max-height:calc(100vh - 120px);overflow-y:auto;padding:6px 10px 8px}
.cra-sec{margin:4px 0 8px}
.cra-h{font-weight:600;text-transform:uppercase;letter-spacing:.04em;font-size:11px;color:#666;margin:6px 0 4px}
.cra-row{margin:3px 0;display:flex;gap:10px;align-items:center;flex-wrap:wrap}
.cra-row select{max-width:200px;font-size:12px}
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
  font:12px/1.3 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;width:250px}
.cra-legend .t{font-weight:600}
.cra-legend .s{color:#555;margin:2px 0 5px}
.cra-legend .bar{height:12px;border-radius:2px}
.cra-legend .ax{display:flex;justify-content:space-between;margin-top:2px}
.cra-legend .n{color:#777;font-size:11px;margin-top:3px}
]---"

js <- r"---[
function(el, x, cfg) {
  var map = this, lm = map.layerManager, idx = {};
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
    var a = Math.abs(v);
    return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) :
      String(Number(v.toPrecision(2)));
  };
  var st = {id: cfg.def_id, type: 'value', period: cfg.def_period, op: cfg.opacity, cur: null};

  A(cfg.layers).forEach(function(g) { on(g, false); });

  var h = '<div class="cra-hd"><span>' + esc(cfg.title) +
    '</span><button type="button" class="cra-tg" title="Collapse">&#8722;</button></div><div class="cra-bd">';
  h += '<div class="cra-sec"><div class="cra-h">Indicator layer</div>' +
    '<div class="cra-row"><label><input type="radio" name="cra-ty" value="value" checked> Values</label>' +
    '<label><input type="radio" name="cra-ty" value="score"> Scores (0&ndash;1)</label></div>' +
    '<div class="cra-row"><label>Period <select class="cra-per">' +
    A(cfg.periods).map(function(p) {
      return '<option value="' + esc(p) + '"' + (p === st.period ? ' selected' : '') + '>' +
        esc(cfg.plabs[p] || p) + '</option>';
    }).join('') + '</select></label></div>' +
    '<div class="cra-row"><label>Opacity <input type="range" class="cra-op" min="0" max="1" step="0.05" value="' +
    st.op + '"></label></div>';
  h += '<label class="cra-ind" style="margin-left:0"><input type="radio" name="cra-ind" value=""> None</label>';
  A(cfg.tree).forEach(function(sr) {
    h += '<div class="cra-sr" title="' + esc(sr.desc) + '">' + esc(sr.code) +
      (sr.desc ? ': ' + esc(sr.desc) : '') + '</div>';
    A(sr.comps).forEach(function(c) {
      h += '<details class="cra-c" open><summary>' + esc(c.name) + '</summary>';
      A(c.factors).forEach(function(f) {
        var fi = A(f.inds), op = fi.some(function(i) { return i.id === st.id; });
        h += '<details class="cra-f"' + (op ? ' open' : '') + '><summary>' + esc(f.name) + '</summary>';
        fi.forEach(function(i) {
          if (i.built) idx[i.id] = i;
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
  h += '</div></div>';

  var panel = L.control({position: 'topright'});
  panel.onAdd = function() {
    var d = L.DomUtil.create('div', 'cra-panel');
    d.innerHTML = h;
    L.DomEvent.disableClickPropagation(d);
    L.DomEvent.disableScrollPropagation(d);
    return d;
  };
  panel.addTo(map);
  var pd = panel.getContainer(), perSel = pd.querySelector('.cra-per');

  var legend = L.control({position: 'bottomleft'});
  legend.onAdd = function() { return L.DomUtil.create('div', 'cra-legend'); };
  legend.addTo(map);
  var ld = legend.getContainer();

  var front = function() { if (st.cur) G(st.cur).eachLayer(function(l) { if (l.bringToFront) l.bringToFront(); }); };
  var setOp = function() { if (st.cur) G(st.cur).eachLayer(function(l) { if (l.setOpacity) l.setOpacity(st.op); }); };

  function draw() {
    if (st.cur) on(st.cur, false);
    st.cur = null;
    var i = idx[st.id];
    if (!i) { ld.style.display = 'none'; perSel.disabled = true; return; }
    var ps = A(i.periods), p = ps.indexOf(st.period) >= 0 ? st.period : ps[0];
    perSel.disabled = ps.length < 2;
    st.cur = [st.id, st.type, p].join('|');
    on(st.cur, true);
    setOp();
    front();
    var lg = cfg.legends[st.id + '|' + st.type];
    var sub = st.type === 'score' ? 'Score: 0 = lowest, 1 = highest risk' :
      'Value' + (i.units ? ' (' + i.units + ')' : '');
    var note = st.type === 'value' ? (i.dir == -1 ? 'Lower values = higher risk' : 'Higher values = higher risk') : '';
    ld.style.display = '';
    ld.innerHTML = '<div class="t">' + esc(i.label) + '</div><div class="s">' + esc(sub) + '<br>' +
      esc(cfg.plabs[p] || p) + '</div><div class="bar" style="background:linear-gradient(to right,' +
      A(lg.cols).join(',') + ')"></div><div class="ax"><span>' + fmt(lg.min) + '</span><span>' +
      fmt(lg.max) + '</span></div>' + (note ? '<div class="n">' + note + '</div>' : '');
  }

  each('input[name=cra-ind]', function(e) { e.addEventListener('change', function() { st.id = this.value; draw(); }); });
  each('input[name=cra-ty]', function(e) { e.addEventListener('change', function() { st.type = this.value; draw(); }); });
  perSel.addEventListener('change', function() { st.period = this.value; draw(); });
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
