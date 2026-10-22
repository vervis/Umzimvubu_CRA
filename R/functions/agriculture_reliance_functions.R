# ============================================================================
# R/functions/agriculture_reliance_functions.R
#
# Vulnerability indicator 203: reliance on agriculture for livelihoods.
#
# Source: Stats SA Census 2011 (full count), agricultural household question
# H14 ("What kind of agricultural activity is this household involved in?"),
# re-tabulated on 2016 wards: SuperWEB2 > Census Boundary Alignment Profiles >
# Municipal Boundaries 2016 > Census 2011 > Agricultural Households
# (dataset id Householdagriculture16_wd), household adjusted.
#
# Input table (exported by hand, 30 Sep 2026):
#   data/raw/census2011/Census2011_Ward2016_ActivityUnion_HH.csv
#   Rows = 2016 wards (EC + KZN). Columns = the six activity tick-boxes
#   crossed (Livestock x Poultry x Vegetable x Other crops x Fodder x Other),
#   each with its tick value and "Not applicable" (= not ticked), i.e. 64
#   combinations (+ a Livestock "Total" block, which is ignored).
#
# From that one table, per ward:
#   households          sum of the 64 combinations
#   share_any_activity  1 - (all six "Not applicable") / households  <- indicator
#   share_<activity>    households ticking that activity / households
#   share_any_crop      vegetables or other crops
#
# Checks done on 30 Sep 2026: households with any activity match the Census
# 2011 agricultural households 10% sample (DataFirst 485) almost exactly
# (Umzimvubu 26,735 vs 26,714; Mbizana 30,532 vs 30,533; Greater Kokstad
# 4,680 vs 4,680).
#
# Denominator: the table's "household adjusted" totals are 2-16% above the
# Census 2011 household counts (e.g. Mbizana 56,179 vs 48,396; Umzimkhulu
# 43,615 vs 42,891), for reasons not documented. Numerators match the
# census, so the table total inflates the denominator unevenly between
# municipalities. denominator = "census2011_sal" (default) instead uses
# Census 2011 SAL household counts (full count; source-of-water.dta), spread
# on SANLC settlement and summed by 2016 ward, as in 202. "table" uses the
# table total.
#
# Needs: sf, terra, haven. Uses SANLC_SETTLEMENT_CLASSES and
# sal_household_counts() (water_access_functions.R) and
# read_census2011_layer() and spread_households() (crop_reliance_functions.R).
# tar_source("R/functions") loads all of them.
# ============================================================================


# Short names for the H14 tick values
H14_ACTIVITY_NAMES <- c(
  "Livestock"      = "livestock",
  "Poultry"        = "poultry",
  "Vegetable"      = "vegetables",
  "Other crops"    = "other_crops",
  "Fodder grazing" = "fodder",
  "Other"          = "other"
)


# ----------------------------------------------------------------------------
# read_superweb_table(): one SuperWEB2 export (csv or xlsx) -> data frame
# ----------------------------------------------------------------------------
# Returns geo (row label), ward_id (8-digit code, NA otherwise), lm_code (NA
# otherwise) and one numeric column per table column. For a crossed table the
# labels of each header row are forward-filled and joined with " | ", in the
# order the fields were nested. Attributes: "geography", "counting".

read_superweb_table <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    sep <- if (sum(grepl(";", lines)) > sum(grepl(",", lines))) ";" else ","
    nf  <- max(utils::count.fields(textConnection(lines), sep = sep, quote = "\""),
               na.rm = TRUE)
    raw <- utils::read.table(text = lines, sep = sep, quote = "\"", header = FALSE,
                             colClasses = "character", fill = TRUE,
                             col.names = paste0("V", seq_len(nf)),
                             comment.char = "", blank.lines.skip = FALSE)
  } else {
    raw <- as.data.frame(readxl::read_excel(path, sheet = 1, col_names = FALSE,
                                            col_types = "text", .name_repair = "minimal"))
  }
  raw <- as.data.frame(lapply(raw, function(x) trimws(as.character(x))),
                       stringsAsFactors = FALSE)
  raw[raw == ""] <- NA

  # "Geography 2016" row: column A (csv) or B (xlsx)
  first2 <- apply(raw[, 1:2, drop = FALSE], 1, function(r) paste(stats::na.omit(r), collapse = ""))
  geo_row <- which(grepl("^Geography (2011|2016)$", first2))[1]
  if (is.na(geo_row)) stop(basename(path), ": no 'Geography 2011/2016' row found.")
  geo_col <- which(grepl("^Geography", unlist(raw[geo_row, 1:2])))[1]
  geography <- raw[geo_row, geo_col]

  data_cols <- (geo_col + 1):ncol(raw)
  sum_row   <- which(grepl("^Default Summation", raw[[1]]))[1]
  head_rows <- (sum_row + 1):(geo_row - 1)
  head_rows <- head_rows[rowSums(!is.na(raw[head_rows, data_cols, drop = FALSE])) > 0]

  lab <- lapply(head_rows, function(i) {
    v <- unlist(raw[i, data_cols])
    for (k in seq_along(v)[-1]) if (is.na(v[k])) v[k] <- v[k - 1]
    v
  })
  labels <- do.call(paste, c(lab, sep = " | "))
  keep <- !is.na(unlist(raw[head_rows[length(head_rows)], data_cols]))
  data_cols <- data_cols[keep]
  labels <- labels[keep]

  body <- raw[(geo_row + 1):nrow(raw), , drop = FALSE]
  body <- body[!is.na(body[[geo_col]]) &
                 !grepl("Statistics South Africa|^Total$", body[[geo_col]]), ]
  num <- function(x) {
    x <- gsub(intToUtf8(160), "", x, fixed = TRUE)   # non-breaking space
    x <- gsub(" ", "", x, fixed = TRUE)
    as.numeric(sub(",", ".", x, fixed = TRUE))
  }
  out <- data.frame(geo = body[[geo_col]], stringsAsFactors = FALSE)
  labs_u <- make.unique(labels)
  for (j in seq_along(data_cols)) out[[labs_u[j]]] <- num(body[[data_cols[j]]])
  out <- out[!duplicated(out$geo), ]

  out$ward_id <- ifelse(grepl("^\\d{8}\\s*:", out$geo), sub("^(\\d{8}).*$", "\\1", out$geo), NA)
  out$lm_code <- ifelse(grepl("^[A-Z]{2,3}\\d{3}\\s*:", out$geo),
                        sub("^([A-Z]{2,3}\\d{3}).*$", "\\1", out$geo), NA)
  attr(out, "geography") <- geography
  attr(out, "counting")  <- raw[which(grepl("^Counting", raw[[1]]))[1], 1]
  out
}


# ----------------------------------------------------------------------------
# activity_union_counts(): households, "none", and each activity, per row
# ----------------------------------------------------------------------------
# tab: read_superweb_table() of the crossed table. Columns whose labels
# include "Total" or "Unspecified" are dropped; each remaining column is one
# combination of ticked / "Not applicable" across the fields.

activity_union_counts <- function(tab, activity_names = H14_ACTIVITY_NAMES) {
  cols  <- setdiff(names(tab), c("geo", "ward_id", "lm_code"))
  parts <- strsplit(sub("\\.\\d+$", "", cols), " | ", fixed = TRUE)
  drop  <- vapply(parts, function(p) any(grepl("^(Total|Unspecified)", p)), logical(1))
  cols  <- cols[!drop]; parts <- parts[!drop]
  nfield <- unique(lengths(parts))
  if (length(nfield) != 1) stop("Inconsistent header depth in the crossed table.")
  if (length(cols) != 2^nfield) {
    stop("Expected ", 2^nfield, " combinations for ", nfield, " fields, found ", length(cols),
         ". Tick only the activity value and 'Not applicable' in each field.")
  }
  ticked <- do.call(rbind, lapply(parts, function(p) !grepl("^Not applicable", p)))
  m <- as.matrix(tab[cols])
  m[is.na(m)] <- 0

  out <- data.frame(geo = tab$geo, ward_id = tab$ward_id, lm_code = tab$lm_code,
                    households_table = rowSums(m),
                    n_any_activity = rowSums(m[, rowSums(ticked) > 0, drop = FALSE]))
  for (k in seq_len(nfield)) {
    lab <- unique(vapply(parts, `[`, "", k)[ticked[, k]])
    nm  <- if (lab %in% names(activity_names)) activity_names[[lab]] else make.names(tolower(lab))
    out[[paste0("n_", nm)]] <- rowSums(m[, ticked[, k], drop = FALSE])
  }
  crop <- vapply(parts, function(p) any(p %in% c("Vegetable", "Other crops")), logical(1))
  out$n_any_crop <- rowSums(m[, crop, drop = FALSE])
  out
}


# ----------------------------------------------------------------------------
# prepare_agriculture_reliance()
# ----------------------------------------------------------------------------
#   union_file     : crossed activity table (csv or xlsx; see header)
#   wards_gpkg     : 2016 ward boundaries (106 output)
#   sanlc_file     : SANLC 2022 raster (CRA extract)
#   template       : the `template_grid` target (1 km, SA Albers)
#   out_gpkg       : ward table with counts and shares
#   out_tif        : indicator (share_any_activity) on the 1 km grid
#   out_tif_shares : optional multi-layer 1 km file, one layer per share
#   denominator    : "census2011_sal" (default) or "table" (see header)
#   mapdata_zip, households_dta : Census 2011 SAL boundaries and a full-count
#                    SAL household table (for denominator = "census2011_sal")
#   lm_check_zip   : optional DataFirst 485 zip (Census 2011 agricultural
#                    households 10% sample) for a municipal check
#   fine           : sub-cells per template cell side (10 -> 100 m)
#
# Ward shares go onto the 1 km grid weighted by SANLC 2022 settlement on the
# 100 m sub-grid, so each 1 km value mixes the wards whose settlements fall in
# it. Cells with no mapped settlement are NA.

prepare_agriculture_reliance <- function(union_file, wards_gpkg, sanlc_file,
                                         template, out_gpkg, out_tif,
                                         out_tif_shares = NULL,
                                         denominator = c("census2011_sal", "table"),
                                         mapdata_zip = NULL,
                                         households_dta = NULL,
                                         lm_check_zip = NULL, fine = 10,
                                         settlement_classes = SANLC_SETTLEMENT_CLASSES) {
  denominator <- match.arg(denominator)

  # 1. Counts from the crossed table ----------------------------------------
  tab <- read_superweb_table(union_file)
  if (attr(tab, "geography") != "Geography 2016") {
    stop(basename(union_file), " is on ", attr(tab, "geography"), "; 2016 wards expected.")
  }
  if (!grepl("household", attr(tab, "counting"), ignore.case = TRUE)) {
    stop(basename(union_file), " is not counted by households: ", attr(tab, "counting"))
  }
  cnt <- activity_union_counts(tab)
  n_cols <- grep("^n_", names(cnt), value = TRUE)
  lm_cnt <- cnt[!is.na(cnt$lm_code), ]
  wd <- cnt[!is.na(cnt$ward_id), c("ward_id", "households_table", n_cols)]

  # 2. Wards and the settlement sub-grid ----------------------------------------
  tmpl  <- terra::rast(template)
  sub   <- terra::disagg(terra::rast(tmpl), fact = fine)
  sanlc <- terra::rast(sanlc_file)
  settle <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
  settle <- terra::project(settle, sub, method = "average")
  settle <- terra::ifel(is.na(settle), 0, settle)

  wards <- sf::st_read(wards_gpkg, quiet = TRUE)
  wards <- sf::st_transform(wards, terra::crs(tmpl))
  wards <- wards[, c("ward_id", "ward_no", "lm_code", "lm_name", "dm_code")]
  wards <- merge(wards, wd, by = "ward_id", all.x = TRUE)
  if (any(is.na(wards$households_table))) {
    warning(sum(is.na(wards$households_table)), " wards in the ROI have no census values.")
  }
  wv <- terra::vect(wards)
  wv$zone <- seq_len(nrow(wv))
  wzone <- terra::rasterize(wv, sub, field = "zone")

  # 3. Denominator ------------------------------------------------------------
  if (denominator == "census2011_sal") {
    if (is.null(mapdata_zip) || is.null(households_dta)) {
      stop("denominator = 'census2011_sal' needs mapdata_zip and households_dta.")
    }
    # Sub-grid over the full extent of the wards (wider than the template),
    # so wards cut by the grid edge still get all their households. Outside
    # the SANLC extract there is no settlement layer; SALs there are spread
    # evenly, which still keeps every household in its own SAL and ward.
    res_sub <- terra::res(sub)
    e <- terra::ext(wv)
    e <- terra::ext(floor(e$xmin / res_sub[1]) * res_sub[1], ceiling(e$xmax / res_sub[1]) * res_sub[1],
                    floor(e$ymin / res_sub[2]) * res_sub[2], ceiling(e$ymax / res_sub[2]) * res_sub[2])
    sub_w <- terra::rast(e, resolution = res_sub, crs = terra::crs(tmpl))
    settle_w <- terra::ifel(is.na(sanlc), NA, sanlc %in% settlement_classes) * 1
    settle_w <- terra::project(settle_w, sub_w, method = "average")
    wzone_w <- terra::rasterize(wv, sub_w, field = "zone")

    ward_bbox <- sf::st_bbox(sf::st_transform(wards, 4326))
    hh <- sal_household_counts(households_dta, "wtrsrs")
    hh$hh_total <- rowSums(hh[grep("^wtrsrs_", names(hh))])
    sal <- read_census2011_layer(mapdata_zip, "SAL_SA_2013", ward_bbox)
    sal$sal_code <- sprintf("%.0f", sal$SAL_CODE)
    sal <- merge(sal[, "sal_code"], hh[c("sal_code", "hh_total")], by = "sal_code", all.x = TRUE)
    sv <- terra::project(terra::vect(sal), terra::crs(tmpl))
    hh_sub <- spread_households(sv, sub_w, settle_w, "hh_total")$surface
    z <- terra::zonal(hh_sub, wzone_w, fun = "sum", na.rm = TRUE)
    wards$households_sal2011 <- NA_real_
    wards$households_sal2011[z[[1]]] <- z[[2]]
    den <- wards$households_sal2011
  } else {
    den <- wards$households_table
  }
  wards$households <- den

  # 4. Shares -------------------------------------------------------------------
  share <- function(n) ifelse(den > 0, pmin(n / den, 1), NA_real_)
  for (nc in n_cols) wards[[sub("^n_", "share_", nc)]] <- share(sf::st_drop_geometry(wards)[[nc]])
  wards$share_capped <- sf::st_drop_geometry(wards)$n_any_activity > den
  wards$area_km2 <- as.numeric(sf::st_area(wards)) / 1e6

  dir.create(dirname(out_gpkg), recursive = TRUE, showWarnings = FALSE)
  sf::st_write(wards, out_gpkg, layer = "agriculture_reliance_wards2016",
               delete_dsn = TRUE, quiet = TRUE)

  # 5. Onto the 1 km grid, weighted by settlement -----------------------------
  wv <- terra::vect(wards)
  to_grid <- function(field) {
    s   <- terra::rasterize(wv, sub, field = field)
    num <- terra::aggregate(terra::ifel(is.na(s), 0, s * settle), fact = fine, fun = "sum")
    den <- terra::aggregate(terra::ifel(is.na(s), 0, settle), fact = fine, fun = "sum")
    r <- terra::mask(terra::ifel(den > 0, num / den, NA), tmpl)
    names(r) <- field
    r
  }
  gdal <- c("COMPRESS=DEFLATE", "TILED=YES")
  r <- to_grid("share_any_activity")
  dir.create(dirname(out_tif), recursive = TRUE, showWarnings = FALSE)
  terra::writeRaster(r, out_tif, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
  outs <- c(out_tif, out_gpkg)

  if (!is.null(out_tif_shares)) {
    fields <- grep("^share_", names(wards), value = TRUE)
    fields <- setdiff(fields, "share_capped")
    rs <- terra::rast(lapply(fields, to_grid))
    terra::writeRaster(rs, out_tif_shares, overwrite = TRUE, datatype = "FLT4S", gdal = gdal)
    outs <- c(outs, out_tif_shares)
  }

  # QA -------------------------------------------------------------------------
  w <- sf::st_drop_geometry(wards)
  message("Denominator: ", denominator,
          " | wards in ROI: ", nrow(w), " | share_any_activity range: ",
          paste(round(range(w$share_any_activity, na.rm = TRUE), 3), collapse = " - "),
          " | wards capped at 1: ", sum(w$share_capped, na.rm = TRUE))
  if (denominator == "census2011_sal") {
    message("Table households / SAL 2011 households (ward median, p10-p90): ",
            round(stats::median(w$households_table / w$households_sal2011, na.rm = TRUE), 3), ", ",
            paste(round(stats::quantile(w$households_table / w$households_sal2011,
                                        c(0.1, 0.9), na.rm = TRUE), 3), collapse = "-"))
  }
  message("Study-area cells with a value: ",
          round(100 * terra::global(!is.na(r), "sum")[[1]] /
                  terra::global(!is.na(tmpl), "sum")[[1]], 1), "%")

  if (!is.null(lm_check_zip) && file.exists(lm_check_zip)) {
    exdir <- file.path(tempdir(), "census2011_agric")
    f <- utils::unzip(lm_check_zip, exdir = exdir)
    f <- f[grepl("\\.dta$", f)][1]
    s <- haven::read_dta(f, col_select = c("H_MUNIC", "AGRICULTURE_10PERCENT_WEIGHT"))
    s$lm_name <- as.character(haven::as_factor(s$H_MUNIC))
    s <- stats::aggregate(AGRICULTURE_10PERCENT_WEIGHT ~ lm_name, data = s, FUN = sum)
    names(s)[2] <- "agric_hh_sample"
    lm_cnt$lm_name <- trimws(sub("^[A-Z]{2,3}\\d{3}\\s*:\\s*", "", lm_cnt$geo))
    chk <- merge(lm_cnt[c("lm_code", "lm_name", "n_any_activity")], s,
                 by.x = tolower("lm_name"), by.y = "lm_name")
    chk <- chk[chk$lm_code %in% unique(w$lm_code), ]
    chk$ratio <- round(chk$n_any_activity / chk$agric_hh_sample, 3)
    message("LM check: households with any activity (table) vs Census 2011 ",
            "agricultural households sample (DataFirst 485); ratio ~1 expected ",
            "where the municipality did not change in 2016:")
    print(chk[order(chk$lm_code), ], row.names = FALSE, digits = 6)
  }
  outs
}
