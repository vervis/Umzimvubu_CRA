# ============================================================================
# R/functions/export_functions.R
#
# Stage-80 export helpers. Only function definitions: tar_source("R/functions")
# loads them.
#
# export_pages(): publishes a saveWidget() map (non-selfcontained) to the
#   GitHub Pages folder (default docs/):
#   1. Copies <name>.html to docs/index.html, rewriting references to the
#      dependency folder <name>_files/ as index_files/.
#   2. Replaces docs/index_files/ with a fresh copy of <name>_files/.
#   3. Ensures docs/.nojekyll exists (so Pages serves the _files folder).
#   Returns all written paths, for a format = "file" target.
#
# export_catalogue(): public indicator catalogue for GitHub Pages, built from
#   CRA_data_sources_Claude.xlsx (sheet Indicators, internal columns dropped),
#   the indicator register and outputs/tables/indicator_summary.csv:
#   - docs/indicators.csv: the public Indicators table (UTF-8, Excel-ready).
#   - docs/indicators.html: self-contained page with search/filters, an
#     overview table, one card per indicator (anchor = Register ID, so the
#     map can link to indicators.html#<id>) grouped Component > Factor, with
#     every column under its xlsx group heading, layer statistics per period,
#     a "View on map" link (index.html#<id>) and the Lookups glossary.
# ============================================================================

export_pages <- function(html, dest = "docs", page = "index") {
  src_lib <- sub("\\.html$", "_files", html)
  if (!file.exists(html) || !dir.exists(src_lib))
    stop("Missing map or its dependency folder: ", html, " / ", src_lib)

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(dest, paste0(page, ".html"))
  lib <- file.path(dest, paste0(page, "_files"))

  txt <- readChar(html, file.size(html), useBytes = TRUE)
  txt <- gsub(paste0(basename(src_lib), "/"), paste0(basename(lib), "/"), txt,
              fixed = TRUE, useBytes = TRUE)
  writeChar(txt, out, eos = NULL, useBytes = TRUE)

  unlink(lib, recursive = TRUE)
  dir.create(lib)
  ok <- file.copy(list.files(src_lib, full.names = TRUE), lib, recursive = TRUE)
  if (!all(ok)) stop("Copying ", src_lib, " to ", lib, " failed")

  nj <- file.path(dest, ".nojekyll")
  if (!file.exists(nj)) file.create(nj)

  c(out, nj, list.files(lib, recursive = TRUE, full.names = TRUE))
}


export_catalogue <- function(xlsx, register, summary, dest = "docs",
                             drop = c("Notes for Claude", "Claude response"),
                             map_page = "index.html",
                             repo = "https://github.com/vervis/Umzimvubu_CRA") {
  rd <- function(...) as.data.frame(readxl::read_excel(xlsx, ..., .name_repair = "minimal"))
  g  <- unlist(rd(sheet = "Indicators", n_max = 1, col_names = FALSE, col_types = "text")[1, ])
  x  <- rd(sheet = "Indicators", skip = 1, col_types = "text")
  g  <- g[seq_len(ncol(x))]
  for (i in seq_along(g)[-1]) if (is.na(g[i])) g[i] <- g[i - 1]
  k  <- !names(x) %in% drop & !is.na(names(x)) & nzchar(names(x))
  x  <- x[!is.na(x[["Register ID"]]), k, drop = FALSE]
  g  <- unname(g[k])
  g[names(x)[seq_along(g)] == "To-do flag"] <- g[match("To-do details", names(x))]

  lk  <- rd(sheet = "Lookups", col_names = FALSE, col_types = "text")
  reg <- utils::read.csv(register, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  sm  <- if (file.exists(summary)) utils::read.csv(summary, stringsAsFactors = FALSE) else
    data.frame(file = character(), layer = character())
  sm  <- sm[grepl("_value\\.tif$", sm$file), ]
  sm$id <- sub("_value\\.tif$", "", sm$file)

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  out_csv  <- file.path(dest, "indicators.csv")
  out_html <- file.path(dest, "indicators.html")
  readr::write_excel_csv(x, out_csv, na = "")

  e  <- function(s) htmltools::htmlEscape(ifelse(is.na(s), "", as.character(s)))
  ln <- function(s) gsub("(https?://[^[:space:]<>\";]*[^[:space:]<>\";.,)])",
                         "<a href=\"\\1\" target=\"_blank\" rel=\"noopener\">\\1</a>", e(s))
  sl <- function(s) gsub("[^a-z0-9]+", "-", tolower(s))
  f3 <- function(v) formatC(signif(v, 3), format = "fg", digits = 3, big.mark = ",")
  pl <- function(p) {
    if (p == "baseline") return("Baseline (static)")
    y <- sub("_", "-", regmatches(p, regexpr("\\d{4}_\\d{4}", p)))
    s <- regmatches(p, regexpr("ssp\\d{3}", p))
    if (!length(s)) return(paste("Historical", y))
    sprintf("%s, SSP%s-%s.%s", y, substr(s, 4, 4), substr(s, 5, 5), substr(s, 6, 6))
  }

  sr  <- lk[grepl("^SR\\d", lk[[1]]), ]
  hdr <- c("Code", "Register ID", "Sub-risk", "Component", "Factor", "Indicator",
           "Question answered", "Status", "Role")
  cmp <- c("Hazard", "Vulnerability", "Exposure")
  x$.c <- match(x$Component, c(cmp, unique(setdiff(x$Component, cmp))))
  x <- x[order(x$.c, match(x$Factor, unique(x$Factor))), ]
  built <- x[["Register ID"]] %in% sm$id
  act   <- x[["Register ID"]] %in% reg$id[as.logical(reg$active)]

  badge <- function(s, cl) sprintf('<span class="b %s">%s</span>', cl, e(s))
  card <- function(i) {
    r  <- x[i, ]
    id <- r[["Register ID"]]
    st <- sm[sm$id == id, ]
    h  <- sprintf(paste0('<article class="card" id="%s" data-c="%s" data-s="%s">',
                         '<div class="ch"><div><div class="fac">%s &rsaquo; %s</div><h3>%s</h3>',
                         '<div class="q">%s</div></div><div class="bs">%s%s%s</div></div>',
                         '<div class="links">%s<a href="#top">Back to overview</a></div>'),
                  e(id), e(r$Component), if (act[i]) "Active" else "Off",
                  e(r$Component), e(r$Factor), e(r$Indicator), e(r[["Question answered"]]),
                  badge(r$Component, sl(r$Component)),
                  badge(r$Status, if (act[i]) "on" else "off"), badge(r$Role, "role"),
                  if (built[i]) sprintf('<a href="%s#%s">View on map &rarr;</a>', map_page, e(id)) else
                    '<span class="na">Not on the map (layer not built)</span>')
    go <- c("Definition", "Rationale & to-do", "Data", "Provenance", "Pipeline")
    for (gg in c(intersect(go, g), setdiff(unique(g), go))) {
      j <- which(g == gg & !names(x)[seq_along(g)] %in% hdr)
      j <- j[!is.na(unlist(r[j])) & nzchar(unlist(r[j]))]
      if (!length(j)) next
      h <- paste0(h, sprintf('<h4>%s</h4><dl>', e(gg)),
                  paste0(sprintf('<dt>%s</dt><dd>%s</dd>', e(names(x)[j]), ln(unlist(r[j]))),
                         collapse = ""), '</dl>')
    }
    h <- paste0(h, sprintf('<dl><dt>Register ID</dt><dd><code>%s</code></dd><dt>Framework code</dt><dd>%s</dd></dl>',
                           e(id), e(r$Code)))
    if (nrow(st)) {
      u <- reg$units[match(id, reg$id)]
      h <- paste0(h, sprintf('<h4>Layer statistics</h4><p class="m">Values%s on the 1 km CRA grid; coverage = share of grid cells with data.</p>',
                             if (!is.na(u)) paste0(" (", e(u), ")") else ""),
                  '<table class="st"><tr><th>Period</th><th>Min</th><th>Mean</th><th>Max</th><th>Coverage</th></tr>',
                  paste0(sprintf('<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s%%</td></tr>',
                                 e(vapply(st$layer, pl, "")), f3(st$min), f3(st$mean), f3(st$max),
                                 formatC(100 * st$coverage, format = "f", digits = 1)),
                         collapse = ""), '</table>')
    }
    paste0(h, '</article>')
  }

  ov <- paste0(sprintf(paste0('<tr data-c="%s" data-s="%s"><td>%s</td><td>%s</td>',
                              '<td><a href="#%s">%s</a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>'),
                       e(x$Component), ifelse(act, "Active", "Off"), e(x$Component), e(x$Factor),
                       e(x[["Register ID"]]), e(x$Indicator), badge(x$Status, ifelse(act, "on", "off")),
                       e(x$Role), e(x$Units), ifelse(built, sprintf('<a href="%s#%s">Map</a>', map_page,
                                                                     e(x[["Register ID"]])), "&ndash;")),
               collapse = "")

  body <- ""
  for (cc in unique(x$Component)) {
    body <- paste0(body, sprintf('<section class="comp" data-c="%s"><h2>%s</h2>', e(cc), e(cc)))
    a <- which(x$Component == cc)
    for (ff in unique(x$Factor[a])) {
      b <- a[x$Factor[a] == ff]
      body <- paste0(body, sprintf('<div class="fg"><h3 class="fh">%s</h3>', e(ff)),
                     paste0(vapply(b, card, ""), collapse = ""), '</div>')
    }
    body <- paste0(body, '</section>')
  }

  gl <- ""
  sec <- NA
  for (i in seq_len(nrow(lk))) {
    a <- lk[i, 1]; b <- if (ncol(lk) > 1) lk[i, 2] else NA
    if (!is.na(a) && is.na(b)) {
      if (!is.na(sec)) gl <- paste0(gl, '</dl>')
      sec <- a
      gl <- paste0(gl, sprintf('<h4>%s</h4><dl>', e(a)))
    } else if (!is.na(b) && !is.na(sec) && !sec %in% c("Notes", "To-do flag")) {
      gl <- paste0(gl, sprintf('<dt>%s</dt><dd>%s</dd>', e(a), ln(b)))
    }
  }
  if (!is.na(sec)) gl <- paste0(gl, '</dl>')

  css <- r"---[
*{box-sizing:border-box}
body{margin:0;font:14px/1.5 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;color:#1d2329;background:#f6f7f9}
a{color:#1f6fb4}
header{background:#1f4e79;color:#fff;padding:22px 16px}
header a{color:#fff}
.w{max-width:1100px;margin:0 auto;padding:0 16px}
header .w{padding:0}
h1{margin:0 0 4px;font-size:24px}
.sub{opacity:.9;max-width:820px}
.nav{margin-top:10px;display:flex;gap:16px;flex-wrap:wrap;font-weight:600}
.ctl{display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:16px 0}
.ctl input[type=search]{flex:1;min-width:200px;padding:7px 10px;border:1px solid #c9cfd6;border-radius:6px;font-size:14px}
.ctl select{padding:6px;border:1px solid #c9cfd6;border-radius:6px}
.tw{overflow-x:auto;background:#fff;border:1px solid #e1e5ea;border-radius:8px}
table{border-collapse:collapse;width:100%}
th,td{text-align:left;padding:6px 10px;border-bottom:1px solid #eef0f3;vertical-align:top}
th{background:#f0f3f6;font-size:12px;text-transform:uppercase;letter-spacing:.03em;color:#4a5560}
h2{margin:28px 0 8px;padding-bottom:4px;border-bottom:2px solid #1f4e79;color:#1f4e79}
.fh{margin:18px 0 6px;font-size:15px;color:#4a5560}
.card{background:#fff;border:1px solid #e1e5ea;border-radius:8px;padding:14px 16px;margin:0 0 12px;scroll-margin-top:12px}
.card:target{border-color:#e66101;box-shadow:0 0 0 3px rgba(230,97,1,.25)}
.ch{display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap}
.card h3{margin:2px 0;font-size:18px}
.fac{font-size:12px;color:#66717c}
.q{font-style:italic;color:#4a5560}
.bs{display:flex;gap:5px;align-items:flex-start;flex-wrap:wrap}
.b{font-size:11px;padding:1px 7px;border-radius:10px;border:1px solid #c9cfd6;white-space:nowrap}
.b.on{background:#e3f4e8;border-color:#7fc596;color:#1d6b35}
.b.off{background:#f0f0f0;color:#666}
.b.hazard{background:#fce4d6;border-color:#f0b48e}
.b.vulnerability{background:#e2efda;border-color:#a9d08e}
.b.exposure{background:#ddebf7;border-color:#9dc3e6}
.links{margin:8px 0 4px;display:flex;gap:16px;font-size:13px;font-weight:600}
.na{color:#8a939c;font-weight:400}
h4{margin:12px 0 4px;font-size:12px;text-transform:uppercase;letter-spacing:.04em;color:#66717c}
dl{display:grid;grid-template-columns:minmax(140px,220px) 1fr;gap:3px 14px;margin:0}
dt{color:#4a5560;font-weight:600}
dd{margin:0;overflow-wrap:anywhere}
.m{color:#66717c;font-size:12px;margin:2px 0 4px}
.st{width:auto;font-size:13px}
.st td,.st th{padding:3px 10px}
.gl{background:#fff;border:1px solid #e1e5ea;border-radius:8px;padding:6px 16px 14px;margin:12px 0 30px}
footer{color:#66717c;font-size:12px;padding:10px 0 30px}
@media (max-width:640px){dl{grid-template-columns:1fr}dt{margin-top:4px}}
]---"

  js <- r"---[
(function(){
  var q=document.getElementById('q'),c=document.getElementById('fc'),s=document.getElementById('fs');
  function run(){
    var t=q.value.toLowerCase(),cv=c.value,sv=s.value;
    var ok=function(el){return (!cv||el.dataset.c===cv)&&(!sv||el.dataset.s===sv)&&(!t||el.textContent.toLowerCase().indexOf(t)>=0);};
    document.querySelectorAll('tbody tr').forEach(function(r){r.style.display=ok(r)?'':'none';});
    document.querySelectorAll('.card').forEach(function(r){r.style.display=ok(r)?'':'none';});
    document.querySelectorAll('.fg').forEach(function(f){var v=[].some.call(f.querySelectorAll('.card'),function(k){return k.style.display!=='none';});f.style.display=v?'':'none';});
    document.querySelectorAll('.comp').forEach(function(f){var v=[].some.call(f.querySelectorAll('.card'),function(k){return k.style.display!=='none';});f.style.display=v?'':'none';});
  }
  [q,c,s].forEach(function(el){el.addEventListener('input',run);el.addEventListener('change',run);});
})();
]---"

  opts <- paste0(sprintf('<option value="%s">%s</option>', e(unique(x$Component)), e(unique(x$Component))),
                 collapse = "")
  html <- paste0(
    '<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>Umzimvubu CRA indicators</title><style>', css, '</style></head><body>',
    '<header id="top"><div class="w"><h1>Umzimvubu CRA: indicator catalogue</h1><div class="sub">',
    if (nrow(sr)) paste0(e(sr[1, 1]), ": ", e(sr[1, 2]), ". ") else "",
    'Details of every indicator layer in the Climate Risk Assessment: definition, data source, ',
    'processing, licence and current status.</div><div class="nav">',
    sprintf('<a href="%s">Open the map</a>', map_page),
    '<a href="indicators.csv" download>Download table (CSV)</a>',
    sprintf('<a href="%s">Code and methods (GitHub)</a>', repo),
    '</div></div></header><main class="w">',
    '<div class="ctl"><input type="search" id="q" placeholder="Search indicators, sources, notes...">',
    '<select id="fc"><option value="">All components</option>', opts, '</select>',
    '<select id="fs"><option value="">All statuses</option><option value="Active">Active</option>',
    '<option value="Off">Off</option></select></div>',
    '<div class="tw"><table><thead><tr><th>Component</th><th>Factor</th><th>Indicator</th>',
    '<th>Status</th><th>Role</th><th>Units</th><th>Map</th></tr></thead><tbody>', ov,
    '</tbody></table></div>', body,
    '<h2>Glossary</h2><div class="gl">', gl, '</div>',
    sprintf('<footer>Generated %s by the targets pipeline (R/functions/export_functions.R) from ',
            format(Sys.Date(), "%d %b %Y")),
    'CRA_data_sources_Claude.xlsx, R/30_indicators/indicator_register.csv and ',
    'outputs/tables/indicator_summary.csv. Status reflects the pipeline at that date. ',
    'Draft research output: indicative only, not for site-level decisions.</footer>',
    '</main><script>', js, '</script></body></html>'
  )
  con <- file(out_html, open = "w", encoding = "UTF-8")
  writeLines(html, con, useBytes = FALSE)
  close(con)
  c(out_html, out_csv)
}
