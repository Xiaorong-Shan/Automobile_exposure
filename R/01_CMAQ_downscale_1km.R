# ============================================================
# Auto PM2.5 (with_filter) Annual (2011–2020) + 1km Downscale
# /home/xshan2/HAQ_LAB/Sumaiya/cmaq/cmaq_output/R/data/annual_rds/pm25_sectors/with_filter
# Inputs: annual_rds/pm25_sectors/with_filter (NRD/ONR daily long tables)
# Roadiness: 1km CSV with normalized weights
# Deps: raster, data.table
# Outputs: GeoTIFF + annual tables + QC txt
# ============================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(raster)
  library(data.table)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
years <- 2011:2020

# where you copied the with_filter rds to (Auto folder)
auto_base <- "/scratch/xshan2/R_Code/disperseR/Auto"
nrd_dir <- file.path(auto_base, "NRD")
onr_dir <- file.path(auto_base, "ONR")

# output directory
out_dir <- file.path(auto_base, "AUTO_annual_downscale_1km")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# roadiness (1km)
road_csv <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/roadiness_1km_hw_loc_sherry.csv"
road_weight_col <- "leng.distm2_hw_norm"   # highway
# road_weight_col <- "leng.distm2_loc_norm" # local
tag <- "hw"

# projection (match roadiness)
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000 +units=m +no_defs"

# coarse grid resolution (m)
cmaq_res <- 12000

# sample size for daily QC (optional)
sample_n <- 2e6

# -----------------------------
# HELPERS
# -----------------------------
parse_ym <- function(f) {
  bn <- basename(f)
  m <- regexpr("(20[0-9][0-9])(0[1-9]|1[0-2])", bn, perl = TRUE)
  if (m[1] != -1) {
    s <- regmatches(bn, m)
    return(list(year = as.integer(substr(s, 1, 4)),
                month = as.integer(substr(s, 5, 6))))
  }
  stop("Cannot parse year/month from filename: ", bn)
}

days_in_month <- function(y, m) {
  if (m == 12) {
    d1 <- as.Date(sprintf("%04d-12-01", y))
    d2 <- as.Date(sprintf("%04d-01-01", y + 1))
  } else {
    d1 <- as.Date(sprintf("%04d-%02d-01", y, m))
    d2 <- as.Date(sprintf("%04d-%02d-01", y, m + 1))
  }
  as.integer(d2 - d1)
}

require_var <- function(nc, v) {
  vnames <- names(nc$var)
  if (!(v %in% vnames)) {
    stop("Variable not found: ", v, "\nAvailable vars include:\n  - ",
         paste(vnames, collapse = "\n  - "))
  }
  v
}

# ---- detect bad "days" using sampling across grid ----
# pm3d: [col,row,t] numeric
# ndays: integer
# return: bad day indices (1..ndays)
detect_bad_days <- function(pm3d, ndays,
                            sample_n = 300000,
                            max_cut = 200,
                            p99_cut = 50) {

  nT <- dim(pm3d)[3]
  if (is.null(nT) || nT < 2) return(integer(0))

  # infer steps per day
  steps_per_day <- nT / ndays
  if (abs(steps_per_day - round(steps_per_day)) > 1e-8) {
    warning("Time steps not divisible by ndays. nT=", nT, " ndays=", ndays,
            ". Falling back to chunking by floor(nT/ndays).")
    steps_per_day <- floor(nT / ndays)
    if (steps_per_day < 1) return(integer(0))
  } else {
    steps_per_day <- as.integer(round(steps_per_day))
  }

  bad <- integer(0)

  for (d in seq_len(ndays)) {
    t1 <- (d - 1) * steps_per_day + 1
    t2 <- min(d * steps_per_day, nT)

    # daily mean over that day's time slice (handles hourly by averaging)
    day_mean <- apply(pm3d[, , t1:t2, drop=FALSE], c(1,2), mean, na.rm=TRUE)

    v <- as.vector(day_mean)
    v <- v[is.finite(v)]
    if (length(v) == 0) next
    if (length(v) > sample_n) v <- sample(v, sample_n)

    p99 <- as.numeric(quantile(v, 0.99, na.rm=TRUE))
    mx  <- as.numeric(max(v, na.rm=TRUE))

    if (is.finite(mx) && (mx > max_cut) || (is.finite(p99) && p99 > p99_cut)) {
      bad <- c(bad, d)
    }
  }

  unique(bad)
}

# ---- replace bad days with last good day (locf) ----
# pm3d: [col,row,t]
# bad_days: indices in 1..ndays
# ndays: integer
# returns fixed pm3d
locf_fill_bad_days <- function(pm3d, bad_days, ndays) {
  nT <- dim(pm3d)[3]
  if (length(bad_days) == 0 || is.null(nT) || nT < 2) return(pm3d)

  steps_per_day <- nT / ndays
  if (abs(steps_per_day - round(steps_per_day)) > 1e-8) {
    steps_per_day <- floor(nT / ndays)
    if (steps_per_day < 1) return(pm3d)
  } else {
    steps_per_day <- as.integer(round(steps_per_day))
  }

  bad_days <- sort(unique(bad_days))
  for (d in bad_days) {
    # find last good day before d
    prev <- setdiff(seq_len(d-1), bad_days)
    if (length(prev) == 0) next
    d0 <- tail(prev, 1)

    t1  <- (d - 1) * steps_per_day + 1
    t2  <- min(d * steps_per_day, nT)
    t10 <- (d0 - 1) * steps_per_day + 1
    t20 <- min(d0 * steps_per_day, nT)

    # overwrite the entire day's time-slices with last good day's slices
    # if different lengths (edge cases), recycle to min length
    L  <- t2 - t1 + 1
    L0 <- t20 - t10 + 1
    Lmin <- min(L, L0)

    pm3d[, , t1:(t1+Lmin-1)] <- pm3d[, , t10:(t10+Lmin-1)]
    if (L > Lmin) {
      # if current day has more slices, repeat last slice of d0
      pm3d[, , (t1+Lmin):t2] <- pm3d[, , (t10+Lmin-1)]
    }
  }

  pm3d
}

# ============================================================
# 0) READ ROADINESS -> 1km RASTER TEMPLATE 
# ============================================================
cat("Reading roadiness CSV:\n  ", road_csv, "\n")

road_dt <- fread(road_csv)
cn <- names(road_dt)
xcol <- if ("x" %in% cn) "x" else if ("X" %in% cn) "X" else stop("No x column in roadiness CSV")
ycol <- if ("y" %in% cn) "y" else if ("Y" %in% cn) "Y" else stop("No y column in roadiness CSV")
stopifnot(road_weight_col %in% cn)

road_dt <- road_dt[, .(x = get(xcol), y = get(ycol), w = get(road_weight_col))]
road_dt[!is.finite(w) | is.na(w), w := 0]

road_r <- rasterFromXYZ(as.data.frame(road_dt), crs = CRS(p4s))
names(road_r) <- paste0("road_", tag)

cat("\n[Roadiness raster template]\n")
print(road_r)

fact <- round(cmaq_res / res(road_r)[1])
if (!is.finite(fact) || fact < 1) stop("Bad aggregation factor. Check roadiness resolution.")
cat("\nAggregation factor (coarse/fine): ", fact, "\n")

road_mean_12km <- aggregate(road_r, fact = fact, fun = mean, na.rm = TRUE)
ex_road <- extent(road_r)

# ============================================================
# 1) INDEX CMAQ FILES 
# ============================================================
cat("\nIndexing CMAQ files under:\n  ", post_dir, "\n")
nc_files <- list.files(post_dir, pattern = "hr2day_w_organics.*\\.nc$", full.names = TRUE, recursive = TRUE)
if (length(nc_files) == 0) stop("No hr2day_w_organics netCDF files found under post_dir.")

meta <- rbindlist(lapply(nc_files, function(f) {
  ym <- try(parse_ym(f), silent = TRUE)
  if (inherits(ym, "try-error")) return(NULL)
  data.table(file = f, year = ym$year, month = ym$month, base = basename(f))
}), fill = TRUE)

meta <- meta[year %in% years]
setorder(meta, year, month, base)

cat("\nFiles found per year (before de-dup):\n")
print(meta[, .N, by = year])

meta <- meta[, .SD[1], by = .(year, month)]
setorder(meta, year, month)

cat("\nFiles used per year (after de-dup by year-month):\n")
print(meta[, .N, by = year])

# ============================================================
# 2) MAIN LOOP: YEARS x VARS (WITH blow-up fix)
# ============================================================
# thresholds for blow-up detection (tune if needed)
sample_n <- 300000
max_cut  <- 200
p99_cut  <- 50

for (yy in years) {

  files_y <- meta[year == yy]
  if (nrow(files_y) == 0) {
    cat("\n[SKIP] No files found for year:", yy, "\n")
    next
  }

  qc_year <- file.path(out_dir, sprintf("QC_blowup_fix_%d_%s.txt", yy, tag))
  cat("QC blow-up fix report | year=", yy, " tag=", tag, "\n\n", file=qc_year)

  cat("\n========================================\n")
  cat("YEAR:", yy, "| months:", paste(sprintf("%02d", files_y$month), collapse = ","), "\n")
  cat("========================================\n")

  for (v in vars) {

    cat("\n  ---- Processing variable:", v, "----\n")
    cat("Variable: ", v, "\n", file=qc_year, append=TRUE)

    sum_mat <- NULL
    tot_days <- 0L
    ncol_cmaq <- NULL
    nrow_cmaq <- NULL

    for (i in 1:nrow(files_y)) {
      f  <- files_y$file[i]
      mo <- files_y$month[i]
      nd <- days_in_month(yy, mo)

      cat("    Month", sprintf("%02d", mo), "| days =", nd, "|", basename(f), "\n")

      nc <- nc_open(f)
      vv <- require_var(nc, v)
      pm <- ncvar_get(nc, vv)
      nc_close(nc)

      # dims: [COL, ROW, T] or [COL, ROW, LAY, T]
      if (length(dim(pm)) == 4) pm <- pm[, , 1, ]  # surface layer

      # ---- NEW: blow-up detection + fix on the time dimension ----
      bad_days <- detect_bad_days(pm, ndays = nd,
                                  sample_n = sample_n,
                                  max_cut  = max_cut,
                                  p99_cut  = p99_cut)

      if (length(bad_days)) {
        cat("      >>> bad days detected: ", paste(bad_days, collapse=","), "\n")
        cat(sprintf("Month %02d: %s | bad_days = %s\n",
                    mo, basename(f), paste(bad_days, collapse=",")),
            file=qc_year, append=TRUE)

        pm <- locf_fill_bad_days(pm, bad_days, ndays = nd)
      } else {
        cat(sprintf("Month %02d: %s | bad_days = none\n",
                    mo, basename(f)),
            file=qc_year, append=TRUE)
      }

      # Monthly mean over time dimension (after fixing)
      pm_mean <- apply(pm, c(1, 2), mean, na.rm = TRUE)

      # enforce physical: negatives -> 0 (optional but consistent with your later rule)
      pm_mean[!is.finite(pm_mean)] <- NA_real_
      pm_mean[pm_mean < 0] <- 0

      if (is.null(sum_mat)) {
        sum_mat <- pm_mean * nd
        ncol_cmaq <- dim(pm_mean)[1]
        nrow_cmaq <- dim(pm_mean)[2]
      } else {
        sum_mat <- sum_mat + pm_mean * nd
      }

      tot_days <- tot_days + nd

      rm(pm, pm_mean); gc()
    }

    pm_annual <- sum_mat / tot_days

    # Build 12km annual raster anchored to roadiness xmin/ymin
    r_cmaq <- raster(
      nrows = nrow_cmaq,
      ncols = ncol_cmaq,
      xmn = ex_road@xmin,
      xmx = ex_road@xmin + ncol_cmaq * cmaq_res,
      ymn = ex_road@ymin,
      ymx = ex_road@ymin + nrow_cmaq * cmaq_res,
      crs = CRS(p4s)
    )

    values(r_cmaq) <- as.vector(t(pm_annual))
    names(r_cmaq)  <- paste0(v, "_annual_", yy)

    if (isTRUE(save_coarse_annual)) {
      out_cmaq <- file.path(out_dir, sprintf("CMAQ_%s_annual_%d_anchored.tif", v, yy))
      writeRaster(r_cmaq, out_cmaq, overwrite = TRUE)
      cat("    Saved 12km annual raster:\n      ", out_cmaq, "\n")
    }

    # --- Mean-preserving downscale to 1km ---
    cmaq_on_1km <- resample(r_cmaq, road_r, method = "ngb")
    road_mean_on_1km <- resample(road_mean_12km, road_r, method = "ngb")

    eps <- 1e-12
    den <- road_mean_on_1km
    den[is.na(den) | den < eps] <- NA

    factor_mp <- road_r / den
    factor_mp[is.na(factor_mp) | !is.finite(factor_mp)] <- 1

    down_1km <- cmaq_on_1km * factor_mp
    names(down_1km) <- sprintf("%s_downscaled1km_%s_%d", v, tag, yy)

    out_down <- file.path(out_dir, sprintf("DOWN1KM_%s_%s_%d.tif", v, tag, yy))
    writeRaster(down_1km, out_down, overwrite = TRUE)

    # Quick sanity checks
    m0 <- cellStats(cmaq_on_1km, "mean")
    m1 <- cellStats(down_1km, "mean")
    mx <- cellStats(down_1km, "max")

    cat("    Saved 1km downscaled raster:\n      ", out_down, "\n")
    cat(sprintf("    mean(cmaq_on_1km)=%.6f | mean(down_1km)=%.6f | max(down_1km)=%.4f\n", m0, m1, mx))

    cat(sprintf("Year %d | %s | mean_coarse=%.6f mean_1km=%.6f max_1km=%.6f\n\n",
                yy, v, cellStats(r_cmaq,"mean"), m1, mx),
        file=qc_year, append=TRUE)
  }
}

cat("\nDONE ✅\nOutputs written to:\n  ", out_dir, "\n")
