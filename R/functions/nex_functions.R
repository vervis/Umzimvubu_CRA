#************************************************************
# nex_functions.R
# Shared helpers for the NEX-GDDP-CMIP6 download check (101b) and
# hazard preparation (201).
#
#   nex_inventory()    : parse every .nc in the period folders, keep one
#                        version (default v2.0) and flag files that are not
#                        real NetCDF (server error pages saved as .nc).
#   nex_availability() : valid years per period x model x variable versus
#                        the years in download_jobs.csv.
#   nex_complete()     : models with complete data for every variable.
#   nex_quarantine()   : move invalid files to <raw_dir>/_failed/<period>/
#                        so that cmip6_dl(force = FALSE) requests them again.
#   nex_todo(), nex_same_grid() : missing files and grid check for re-downloads.
#   nex_refetch()      : re-request missing/invalid files from NCSS with
#                        retries (101c).
#   nex_aws_index(), nex_aws_fetch() : fill missing files from the AWS Open
#                        Data copy (101d).
#************************************************************

NEX_RX <- paste0(
  "^([a-z]+)_day_(.+)_(historical|ssp[0-9]{3})_",
  "(r[0-9]+i[0-9]+p[0-9]+f[0-9]+)_([a-z0-9]+)_([0-9]{4})(_v[0-9.]+)?\\.nc$"
)

is_netcdf <- function(f) {
  vapply(f, function(x) {
    b <- readBin(x, "raw", 4)
    length(b) == 4 &&
      (identical(b, as.raw(c(0x89, 0x48, 0x44, 0x46))) ||
         identical(b[1:3], charToRaw("CDF")))
  }, logical(1), USE.NAMES = FALSE)
}

nex_inventory <- function(raw_dir, periods, vars, models, version = "v2.0") {
  f <- unlist(lapply(file.path(raw_dir, periods), list.files,
                     pattern = "\\.nc$", full.names = TRUE))
  bn <- basename(f)
  m <- regmatches(bn, regexec(NEX_RX, bn))
  ok <- lengths(m) == 8
  f <- f[ok]; m <- m[ok]
  inv <- tibble::tibble(
    file     = f,
    period   = basename(dirname(f)),
    variable = vapply(m, `[`, "", 2),
    model    = vapply(m, `[`, "", 3),
    scenario = vapply(m, `[`, "", 4),
    year     = as.integer(vapply(m, `[`, "", 7)),
    version  = sub("^_", "", vapply(m, `[`, "", 8)),
    bytes    = file.size(f)
  )
  inv <- inv[inv$version == version & inv$variable %in% vars &
               inv$model %in% models, ]
  inv$valid <- inv$bytes > 1000 & is_netcdf(inv$file)
  inv
}

nex_availability <- function(inv, jobs, vars, models) {
  grid <- tidyr::expand_grid(
    jobs[, c("period_name", "scenario", "first_year", "last_year")],
    model = models, variable = vars
  )
  ok <- inv[inv$valid, ]
  res <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    yrs <- seq.int(g$first_year, g$last_year)
    have <- unique(ok$year[ok$period == g$period_name & ok$model == g$model &
                             ok$variable == g$variable])
    bad <- unique(inv$year[!inv$valid & inv$period == g$period_name &
                             inv$model == g$model & inv$variable == g$variable])
    miss <- setdiff(yrs, have)
    tibble::tibble(
      n_expected = length(yrs),
      n_valid = sum(yrs %in% have),
      n_invalid = length(intersect(bad, miss)),
      n_absent = length(setdiff(miss, bad)),
      complete = !length(miss),
      missing_years = paste(miss, collapse = ";")
    )
  })
  dplyr::bind_cols(grid, dplyr::bind_rows(res))
}

nex_complete <- function(av, period) {
  x <- av[av$period_name == period, ]
  if (!nrow(x)) return(character())
  names(which(tapply(x$complete, x$model, all)))
}

nex_quarantine <- function(inv, raw_dir) {
  bad <- inv[!inv$valid, ]
  if (!nrow(bad)) return(invisible(bad))
  dest <- file.path(raw_dir, "_failed", bad$period, basename(bad$file))
  for (d in unique(dirname(dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  moved <- file.rename(bad$file, dest)
  if (any(!moved)) warning(sum(!moved), " invalid files could not be moved")
  invisible(bad[moved, ])
}


# ----------------------------------------------------------------------------
# Shared by the re-download routes (101c NCSS, 101d AWS)
#   nex_todo()      : one row per missing year x model x variable in av.
#   nex_same_grid() : file is readable, >= 360 days, same rows, columns and
#                     extent as ref (a valid file already in raw_dir).
# ----------------------------------------------------------------------------

nex_todo <- function(av) {
  miss <- av[!av$complete & nzchar(av$missing_years), ]
  if (!nrow(miss)) return(NULL)
  do.call(rbind, lapply(seq_len(nrow(miss)), function(i) {
    m <- miss[i, ]
    data.frame(period = m$period_name, scenario = m$scenario, model = m$model,
               variable = m$variable,
               year = as.integer(strsplit(m$missing_years, ";")[[1]]))
  }))
}

nex_same_grid <- function(f, ref) {
  r <- terra::rast(ref)
  e <- as.vector(terra::ext(r))
  vapply(f, function(x) {
    g <- tryCatch(terra::rast(x), error = function(err) NULL)
    !is.null(g) && terra::nrow(g) == terra::nrow(r) && terra::ncol(g) == terra::ncol(r) &&
      terra::nlyr(g) >= 360 && max(abs(as.vector(terra::ext(g)) - e)) < 1e-3
  }, logical(1), USE.NAMES = FALSE)
}


# ----------------------------------------------------------------------------
# nex_refetch(): re-request missing/invalid files directly from NCSS
# ----------------------------------------------------------------------------
# The NCCS THREDDS server fails intermittently (HTTP 500, Tomcat
# "IllegalStateException: Failed to obtain an entry from the archive"); the
# same request usually succeeds on a later try. cmip6_dl() saves the error
# page as the .nc and does not retry, so this function retries per file.
#   av       : nex_availability() table (missing_years per period/model/var)
#   ref_file : a valid file whose grid defines the request box
#   rounds   : passes over the files still missing; pause seconds between
#   workers  : parallel requests (httr2 >= 1.1.0), else sequential
#   tries    : immediate retries per request on HTTP >= 500
# Failures look node-specific (some back-end nodes return the error, others
# succeed, ~40-60% success on 21 Oct), so many quick retries work better than
# long pauses.
# Without a time range NCSS returns only the time step nearest to now, so
# each request asks for 1 Jan of the year to 1 Jan of the next (valid in
# noleap, 360-day and standard calendars; each file holds one year). A
# response is kept only if it passes is_netcdf() and nex_same_grid().

NEX_NCSS <- "https://ds.nccs.nasa.gov/thredds/ncss/grid/AMES/NEX/GDDP-CMIP6"

nex_refetch <- function(av, raw_dir, ref_file, rounds = 30, pause = 10,
                        workers = 4, tries = 5, version = "v2.0") {
  todo <- nex_todo(av)
  if (is.null(todo)) return(invisible(NULL))

  known <- basename(list.files(raw_dir, "\\.nc$", recursive = TRUE))
  km <- regmatches(known, regexec(NEX_RX, known))
  km <- do.call(rbind, lapply(km[lengths(km) == 8], function(x)
    data.frame(model = x[3], scenario = x[4], member = x[5], grid = x[6])))
  km <- unique(km)
  pick <- function(mo, sc) {
    k <- km[km$model == mo & km$scenario == sc, ]
    if (!nrow(k)) k <- km[km$model == mo, ]
    if (!nrow(k)) stop("No file name found for model ", mo)
    k[1, c("member", "grid")]
  }
  mg <- unique(todo[c("model", "scenario")])
  mg <- cbind(mg, do.call(rbind, Map(pick, mg$model, mg$scenario)))
  todo <- merge(todo, mg, by = c("model", "scenario"))

  r  <- terra::rast(ref_file)
  e  <- as.vector(terra::ext(r)); h <- terra::res(r) / 2
  bb <- sprintf("&north=%.4f&west=%.4f&east=%.4f&south=%.4f",
                e[[4]] - h[2] + 0.01, e[[1]] + h[1] - 0.01,
                e[[2]] - h[1] + 0.01, e[[3]] + h[2] - 0.01)

  todo$file <- sprintf("%s_day_%s_%s_%s_%s_%d_%s.nc", todo$variable, todo$model,
                       todo$scenario, todo$member, todo$grid, todo$year, version)
  todo$url  <- sprintf(paste0("%s/%s/%s/%s/%s/%s?var=%s%s&horizStride=1",
                              "&time_start=%d-01-01T00:00:00Z&time_end=%d-01-01T00:00:00Z",
                              "&accept=netcdf4-classic"),
                       NEX_NCSS, todo$model, todo$scenario, todo$member, todo$variable,
                       todo$file, todo$variable, bb, todo$year, todo$year + 1L)
  todo$dest <- file.path(raw_dir, todo$period, todo$file)
  for (d in unique(dirname(todo$dest))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  todo$valid <- file.exists(todo$dest)
  todo$valid[todo$valid] <- is_netcdf(todo$dest[todo$valid])

  par <- utils::packageVersion("httr2") >= "1.1.0" && workers > 1
  for (k in seq_len(rounds)) {
    i <- which(!todo$valid)
    if (!length(i)) break
    message("Round ", k, ": requesting ", length(i), " files")
    tmp <- sub("\\.nc$", "_part.nc", todo$dest[i])
    reqs <- lapply(todo$url[i], function(u)
      httr2::request(u) |> httr2::req_timeout(300) |>
        httr2::req_retry(max_tries = tries, backoff = function(i) 1,
                         is_transient = function(resp) httr2::resp_status(resp) >= 500) |>
        httr2::req_error(is_error = function(resp) FALSE))
    if (par) {
      httr2::req_perform_parallel(reqs, paths = tmp, on_error = "continue",
                                  max_active = workers, progress = TRUE)
    } else {
      for (j in seq_along(reqs)) {
        try(httr2::req_perform(reqs[[j]], path = tmp[j]), silent = TRUE)
      }
    }
    ok <- file.exists(tmp) & file.size(tmp) > 1000
    ok[ok] <- is_netcdf(tmp[ok])
    ok[ok] <- nex_same_grid(tmp[ok], ref_file)
    file.rename(tmp[ok], todo$dest[i][ok])
    unlink(tmp[!ok])
    todo$valid[i[ok]] <- TRUE
    message("  valid: ", sum(ok), " | still missing: ", sum(!todo$valid))
    if (k < rounds && any(!todo$valid)) Sys.sleep(pause)
  }
  invisible(todo)
}


# ----------------------------------------------------------------------------
# AWS Open Data route (101d): s3://nex-gddp-cmip6, public HTTPS, no account
# ----------------------------------------------------------------------------
# The S3 files are global (1440 x 600 cells) and chunked one full day per
# chunk, so a spatial subset cannot be read remotely: each file (~180-250 MB)
# is downloaded whole, checked against the md5 in index_v2.0_md5.txt,
# cropped to ref_file's grid, written as a small NetCDF with the same name as
# the NCSS files and deleted.
#   nex_aws_index() : download (once) and parse the v2.0 md5 index.
#   nex_aws_fetch() : todo rows from nex_todo(); workers = parallel downloads;
#                     tmp_dir needs ~250 MB per worker.

NEX_S3 <- "https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com"

nex_aws_index <- function(man_dir, refresh = FALSE) {
  f <- file.path(man_dir, "aws_index_v2.0_md5.txt")
  if (refresh || !file.exists(f))
    httr2::request(paste0(NEX_S3, "/index_v2.0_md5.txt")) |>
      httr2::req_retry(max_tries = 5) |> httr2::req_perform(path = f)
  x <- utils::read.table(f, col.names = c("md5", "key"), colClasses = "character")
  bn <- basename(x$key)
  m <- regmatches(bn, regexec(NEX_RX, bn))
  ok <- lengths(m) == 8
  x <- x[ok, ]; m <- m[ok]
  data.frame(md5 = x$md5, key = x$key, file = bn[ok],
             variable = vapply(m, `[`, "", 2), model = vapply(m, `[`, "", 3),
             scenario = vapply(m, `[`, "", 4),
             year = as.integer(vapply(m, `[`, "", 7)))
}

nex_aws_crop <- function(src, dest, var, ref) {
  r <- terra::rast(src)
  if (terra::xmax(r) > 180 && terra::xmin(ref) < 0) r <- terra::rotate(r)
  x <- terra::crop(r, terra::ext(ref), snap = "near")
  u <- terra::units(r)[1]
  terra::writeCDF(x, dest, varname = var, longname = terra::longnames(r)[1],
                  unit = if (is.na(u)) "" else u, zname = "time",
                  compression = 4, overwrite = TRUE)
  dest
}

nex_aws_fetch <- function(todo, idx, raw_dir, ref_file, workers = 2,
                          tmp_dir = tempdir(), timeout = 3600) {
  todo <- merge(todo, idx, by = c("model", "scenario", "variable", "year"),
                all.x = TRUE)
  todo$dest <- file.path(raw_dir, todo$period, todo$file)
  todo$status <- ifelse(is.na(todo$key), "not on S3", "pending")
  todo$mb <- NA_real_
  ref <- terra::rast(ref_file)
  for (d in unique(dirname(todo$dest[!is.na(todo$key)])))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  par <- utils::packageVersion("httr2") >= "1.1.0" && workers > 1
  i_all <- which(todo$status == "pending")
  t0 <- Sys.time()
  for (b in split(i_all, ceiling(seq_along(i_all) / max(1, workers)))) {
    tmp <- file.path(tmp_dir, todo$file[b])
    reqs <- lapply(todo$key[b], function(k)
      httr2::request(paste0(NEX_S3, "/", k)) |> httr2::req_timeout(timeout) |>
        httr2::req_retry(max_tries = 4))
    if (par) {
      httr2::req_perform_parallel(reqs, paths = tmp, on_error = "continue",
                                  max_active = workers)
    } else {
      for (j in seq_along(reqs))
        try(httr2::req_perform(reqs[[j]], path = tmp[j]), silent = TRUE)
    }
    for (j in seq_along(b)) {
      k <- b[j]
      if (!file.exists(tmp[j])) { todo$status[k] <- "download failed"; next }
      todo$mb[k] <- file.size(tmp[j]) / 2^20
      if (unname(tools::md5sum(tmp[j])) != todo$md5[k]) {
        todo$status[k] <- "md5 mismatch"
      } else {
        part <- sub("\\.nc$", "_part.nc", todo$dest[k])
        res <- tryCatch(nex_aws_crop(tmp[j], part, todo$variable[k], ref),
                        error = function(e) NA_character_)
        if (!is.na(res) && nex_same_grid(part, ref_file)) {
          file.rename(part, todo$dest[k])
          todo$status[k] <- "ok"
        } else {
          unlink(part)
          todo$status[k] <- "crop failed"
        }
      }
      unlink(tmp[j])
    }
    n <- sum(todo$status != "pending" & !is.na(todo$key))
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("%d/%d files | ok %d | %.1f min elapsed, ~%.0f min left",
                    n, length(i_all), sum(todo$status == "ok"), el,
                    el / n * (length(i_all) - n)))
  }
  invisible(todo)
}
