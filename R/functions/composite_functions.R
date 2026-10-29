# ============================================================================
# R/functions/composite_functions.R
#
# Stages 40 (components) and 50 (risk), per sub-risk in the indicator
# register, propagating GCM uncertainty from the per-GCM indicator scores
# written by make_indicator() (stage 30).
#
# make_components(): component index = weighted mean of indicator scores
#   (register `weight`, default 1), ignoring NA indicators; a cell is NA
#   when the available weight is < min_w of the component total.
#   - components with only static indicators: one layer
#     <out>/<subrisk>_<component>_index.tif
#   - components with NEX indicators (hazard): one index per GCM and period
#     <out>/gcm/<subrisk>_<component>_index_<period>.tif (INT2U, x GCM_SCALE)
#     and ensemble statistics across GCMs
#     <out>/<subrisk>_<component>_index_ensemble.tif (<period>_median|p10|p90)
#   - <out>/<subrisk>_components.csv: the index files and their indicators.
#
# make_risk(): risk per GCM and period from the component indices
#   (method "geomean" = (H*V*E)^(1/3), "product" or "mean"), then per period
#   <out>/<subrisk>_risk_ensemble.tif layers
#     <period>_median, <period>_p10, <period>_p90 (across GCMs) and
#     <period>_high_share (% of GCMs with risk >= threshold), where the
#     threshold is the q_high quantile of the historical median risk over the
#     CRA extent (fixed, so shares show change and spatial pattern);
#   <out>/<subrisk>_risk_settings.csv: method, q_high, threshold, periods.
# ============================================================================

COMPONENTS <- c("hazard", "vulnerability", "exposure")


row_q <- function(m, probs) {
  s <- if (ncol(m) > 1) t(apply(m, 1, sort, na.last = TRUE)) else m
  n <- rowSums(!is.na(m))
  i <- which(n > 0)
  vapply(probs, function(p) {
    h <- (n[i] - 1) * p + 1
    lo <- floor(h)
    out <- rep(NA_real_, nrow(m))
    out[i] <- s[cbind(i, lo)] + (h - lo) * (s[cbind(i, pmin(lo + 1, n[i]))] - s[cbind(i, lo)])
    out
  }, numeric(nrow(m)))
}

ens_stats <- function(r, probs = c(0.5, 0.1, 0.9)) {
  out <- terra::rast(r, nlyrs = length(probs))
  terra::values(out) <- row_q(terra::values(r), probs)
  out
}

wmean <- function(rr, w, min_w) {
  num <- 0
  den <- 0
  for (k in seq_along(rr)) {
    ok <- !is.na(rr[[k]])
    num <- num + w[k] * terra::ifel(ok, rr[[k]], 0)
    den <- den + w[k] * ok
  }
  terra::ifel(den >= min_w * sum(w), num / den, NA)
}


make_components <- function(register, indicator_files, periods, models,
                            out_dir = "data/processed/components",
                            ind_dir = "data/processed/indicators", min_w = 0.5) {
  force(indicator_files)
  reg <- read_indicator_register(register, syms = FALSE)
  reg <- reg[reg$component %in% COMPONENTS & reg$weight > 0, ]
  dir.create(file.path(out_dir, "gcm"), recursive = TRUE, showWarnings = FALSE)
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  hp <- periods[grepl("^hist", periods)]
  fp <- setdiff(periods, hp)
  pers <- c(hp, fp)

  out <- character()
  for (s in unique(reg$subrisk)) {
    tab <- list()
    for (cp in COMPONENTS) {
      a <- reg[reg$subrisk == s & reg$component == cp, ]
      if (!nrow(a)) stop("Sub-risk ", s, " has no active ", cp, " indicators.")
      sf <- file.path(ind_dir, cp, paste0(a$id, "_score.tif"))
      dyn <- any(a$temporal == "nex")
      if (!dyn) {
        f <- file.path(out_dir, paste0(s, "_", cp, "_index.tif"))
        x <- wmean(lapply(sf, function(z) terra::rast(z)[[1]]), a$weight, min_w)
        names(x) <- "baseline"
        terra::writeRaster(x, f, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
        out <- c(out, f)
        tab[[cp]] <- data.frame(subrisk = s, component = cp, dynamic = FALSE, file = f,
                                indicators = paste(a$id, collapse = ";"),
                                weights = paste(a$weight, collapse = ";"))
        next
      }
      ens <- list()
      for (p in pers) {
        rr <- lapply(seq_len(nrow(a)), function(k) {
          if (a$temporal[k] == "static") return(terra::rast(sf[k])[[1]])
          read_gcm_scores(file.path(ind_dir, cp, "gcm", paste0(a$id[k], "_score_", p, ".tif")))
        })
        g <- terra::rast(lapply(models, function(m) {
          wmean(lapply(rr, function(z) if (terra::nlyr(z) > 1) z[[m]] else z), a$weight, min_w)
        }))
        names(g) <- models
        f <- file.path(out_dir, "gcm", paste0(s, "_", cp, "_index_", p, ".tif"))
        terra::writeRaster(round(g * GCM_SCALE), f, overwrite = TRUE, datatype = "INT2U",
                           gdal = c(gdal, "PREDICTOR=2"))
        out <- c(out, f)
        e <- ens_stats(g)
        names(e) <- paste0(p, c("_median", "_p10", "_p90"))
        ens[[p]] <- e
      }
      f <- file.path(out_dir, paste0(s, "_", cp, "_index_ensemble.tif"))
      terra::writeRaster(terra::rast(unname(ens)), f, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
      out <- c(out, f)
      tab[[cp]] <- data.frame(subrisk = s, component = cp, dynamic = TRUE, file = f,
                              indicators = paste(a$id, collapse = ";"),
                              weights = paste(a$weight, collapse = ";"))
    }
    f <- file.path(out_dir, paste0(s, "_components.csv"))
    tab <- do.call(rbind, tab)
    tab$min_w <- min_w
    utils::write.csv(tab, f, row.names = FALSE)
    out <- c(out, f)
  }
  out
}


make_risk <- function(component_files, periods, models, method = "geomean",
                      q_high = 0.8, out_dir = "data/processed/risk") {
  if (!method %in% c("geomean", "product", "mean")) stop("Unknown risk method: ", method)
  fun <- switch(method,
                geomean = function(h, v, e) (h * v * e)^(1 / 3),
                product = function(h, v, e) h * v * e,
                mean    = function(h, v, e) (h + v + e) / 3)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  hp <- periods[grepl("^hist", periods)]
  pers <- c(hp, setdiff(periods, hp))

  out <- character()
  for (cf in component_files[grepl("_components\\.csv$", component_files)]) {
    tab <- utils::read.csv(cf, stringsAsFactors = FALSE)
    s <- tab$subrisk[1]
    cdir <- dirname(cf)
    comp <- function(cp, p) {
      x <- tab[tab$component == cp, ]
      if (!x$dynamic) return(terra::rast(x$file)[[1]])
      read_gcm_scores(file.path(cdir, "gcm", paste0(s, "_", cp, "_index_", p, ".tif")))
    }
    rk <- lapply(stats::setNames(pers, pers), function(p) {
      h <- comp("hazard", p)
      v <- comp("vulnerability", p)
      e <- comp("exposure", p)
      fun(h, v, e)
    })
    st <- lapply(rk, ens_stats)
    thr <- stats::quantile(terra::values(st[[hp]][[1]], mat = FALSE), q_high,
                           na.rm = TRUE, names = FALSE)
    lay <- lapply(pers, function(p) {
      hs <- 100 * terra::mean(rk[[p]] >= thr)
      x <- c(st[[p]], hs)
      names(x) <- paste0(p, c("_median", "_p10", "_p90", "_high_share"))
      x
    })
    f <- file.path(out_dir, paste0(s, "_risk_ensemble.tif"))
    terra::writeRaster(terra::rast(lay), f, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
    fs <- file.path(out_dir, paste0(s, "_risk_settings.csv"))
    utils::write.csv(data.frame(subrisk = s, method = method, q_high = q_high,
                                threshold = thr, hist_period = hp,
                                n_models = length(models)), fs, row.names = FALSE)
    out <- c(out, f, fs)
  }
  out
}
