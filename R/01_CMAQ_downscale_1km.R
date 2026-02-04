# ============================================================
# Auto PM2.5 (with_filter) Annual (2011–2020) + 1km Downscale
# Inputs: /scratch/xshan2/R_Code/disperseR/Auto/NRD + /ONR
# /home/xshan2/HAQ_LAB/Sumaiya/cmaq/cmaq_output/R/data/annual_rds/pm25_sectors/with_filter
#   - PM25_TOT_NRD_cmaq_YYYY-01_YYYY-12.rds  (daily long table: x,y,Date,val)
#   - PM25_TOT_ONR_cmaq_YYYY-01_YYYY-12.rds
# Roadiness: 1km CSV with normalized weights
# Deps: raster, data.table
# Outputs (per year):
#   - AUTO_TOTAL_annual_YYYY_anchored.tif         (12km)
#   - DOWN1KM_AUTO_TOTAL_<tag>_YYYY.tif           (1km)
#   - TOTAL_YYYY_annual_xy.rds                    (x,y,total,nrd,onr)
#   - QC_AUTO_TOTAL_<tag>_YYYY.txt                (QC report)
# ============================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(raster)
  library(data.table)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
years <- 2011:2020

auto_base <- "/scratch/xshan2/R_Code/disperseR/Auto"
nrd_dir <- file.path(auto_base, "NRD")
onr_dir <- file.path(auto_base, "ONR")

out_dir <- file.path(auto_base, "AUTO_annual_downscale_1km")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

road_csv <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/roadiness_1km_hw_loc_sherry.csv"
road_weight_col <- "leng.distm2_hw_norm"   # highway
# road_weight_col <- "leng.distm2_loc_norm" # local
tag <- "hw"

p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000 +units=m +no_defs"
cmaq_res <- 12000

# daily QC sampling (bigger = more accurate, slower)
sample_n <- 2e6

# blow-up thresholds (very conservative; catches true explosions)
max_cut <- 200
p99_cut <- 50

# netCDF fill value threshold (critical for 2014-like years)
fill_cut <- -1e30

# -----------------------------
# HELPERS
# -----------------------------
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

# sample-based daily QC
daily_stats_sample_dt <- function(dt, valcol, sample_n=2e6) {
  setDT(dt)
  dt[, {
    v <- get(valcol)
    v <- v[is.finite(v)]
    if (length(v) == 0) return(list(p99=NA_real_, mx=NA_real_))
    if (length(v) > sample_n) v <- sample(v, sample_n)
    list(
      p99 = as.numeric(quantile(v, 0.99, na.rm=TRUE)),
      mx  = as.numeric(max(v, na.rm=TRUE))
    )
  }, by = Date]
}

# LOCF fill for bad days (vectorized per cell)
locf_fill_bad_days_dt <- function(dt, valcol, bad_dates) {
  setDT(dt)
  if (length(bad_dates) == 0) return(dt)

  dt[, Date := as.Date(Date)]
  bad_dates <- as.Date(bad_dates)

  # mark bad -> NA, then locf by (x,y)
  dt[Date %in% bad_dates, (valcol) := NA_real_]
  setorder(dt, x, y, Date)
  dt[, (valcol) := nafill(get(valcol), type="locf"), by=.(x,y)]
  dt
}

# build coarse raster from annual table (x,y as index grid; anchored to road extent)
annual_xy_to_raster <- function(tot_dt, ex_road, cmaq_res, p4s) {
  setDT(tot_dt)
  ux <- sort(unique(tot_dt$x))
  uy <- sort(unique(tot_dt$y))
  ncol_cmaq <- length(ux)
  nrow_cmaq <- length(uy)

  mdt <- dcast(tot_dt, y ~ x, value.var = "pm25_total")
  mat <- as.matrix(mdt[, -"y"])

  r <- raster(
    nrows = nrow_cmaq,
    ncols = ncol_cmaq,
    xmn = ex_road@xmin,
    xmx = ex_road@xmin + ncol_cmaq * cmaq_res,
    ymn = ex_road@ymin,
    ymx = ex_road@ymin + nrow_cmaq * cmaq_res,
    crs = CRS(p4s)
  )
  values(r) <- as.vector(t(mat))
  r
}

# -----------------------------
# 0) ROADINESS TEMPLATE
# -----------------------------
cat("Reading roadiness CSV:\n  ", road_csv, "\n")
road_dt <- fread(road_csv)
stopifnot(road_weight_col %in% names(road_dt))

xcol <- if ("x" %in% names(road_dt)) "x" else if ("X" %in% names(road_dt)) "X" else stop("No x column in roadiness CSV")
ycol <- if ("y" %in% names(road_dt)) "y" else if ("Y" %in% names(road_dt)) "Y" else stop("No y column in roadiness CSV")

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

# -----------------------------
# 1) MAIN LOOP (years)
# -----------------------------
for (yy in years) {

  cat("\n========================================\n")
  cat("YEAR:", yy, "\n")
  cat("========================================\n")

  pat <- sprintf("%d-01_%d-12", yy, yy)
  f_nrd <- file.path(nrd_dir, sprintf("PM25_TOT_NRD_cmaq_%s.rds", pat))
  f_onr <- file.path(onr_dir, sprintf("PM25_TOT_ONR_cmaq_%s.rds", pat))
  stopifnot(file.exists(f_nrd), file.exists(f_onr))

  qc_report <- file.path(out_dir, sprintf("QC_AUTO_TOTAL_%s_%d.txt", tag, yy))
  cat(sprintf("QC report | year=%d tag=%s\nNRD=%s\nONR=%s\n\n", yy, tag, f_nrd, f_onr),
      file = qc_report)

  # ---- read daily tables ----
  nrd <- readRDS(f_nrd); setDT(nrd)
  onr <- readRDS(f_onr); setDT(onr)

  # ---- standardize Date & handle fill ----
  nrd[, Date := as.Date(Date)]
  onr[, Date := as.Date(Date)]

  # fill -> NA  (CRITICAL: fixes 2014-like behavior)
  nrd[PM25_TOT_NRD < fill_cut, PM25_TOT_NRD := NA_real_]
  onr[PM25_TOT_ONR < fill_cut, PM25_TOT_ONR := NA_real_]

  # physical: negatives -> 0 (after fill handled)
  nrd[is.finite(PM25_TOT_NRD) & PM25_TOT_NRD < 0, PM25_TOT_NRD := 0]
  onr[is.finite(PM25_TOT_ONR) & PM25_TOT_ONR < 0, PM25_TOT_ONR := 0]

  # ---- daily QC + identify blow-up dates ----
  stats_nrd <- daily_stats_sample_dt(nrd, "PM25_TOT_NRD", sample_n = sample_n)
  stats_onr <- daily_stats_sample_dt(onr, "PM25_TOT_ONR", sample_n = sample_n)

  bad_nrd <- stats_nrd[(mx > max_cut) | (p99 > p99_cut), Date]
  bad_onr <- stats_onr[(mx > max_cut) | (p99 > p99_cut), Date]
  bad_all <- sort(unique(c(bad_nrd, bad_onr)))

  cat("Bad days (NRD): ",
      if (length(bad_nrd)) paste(bad_nrd, collapse=", ") else "(none)",
      "\n", file=qc_report, append=TRUE)
  cat("Bad days (ONR): ",
      if (length(bad_onr)) paste(bad_onr, collapse=", ") else "(none)",
      "\n", file=qc_report, append=TRUE)
  cat("Bad days (union): ",
      if (length(bad_all)) paste(bad_all, collapse=", ") else "(none)",
      "\n\n", file=qc_report, append=TRUE)

  # ---- LOCF fill bad days (use union so NRD/ONR aligned in time) ----
  if (length(bad_all)) {
    nrd <- locf_fill_bad_days_dt(nrd, "PM25_TOT_NRD", bad_all)
    onr <- locf_fill_bad_days_dt(onr, "PM25_TOT_ONR", bad_all)
  }

  # ---- annual mean by cell ----
  pm_nrd <- nrd[, .(pm25_nrd = mean(PM25_TOT_NRD, na.rm=TRUE)), by=.(x,y)]
  pm_onr <- onr[, .(pm25_onr = mean(PM25_TOT_ONR, na.rm=TRUE)), by=.(x,y)]
  setkey(pm_nrd, x,y); setkey(pm_onr, x,y)

  tot <- pm_nrd[pm_onr]
  tot[, pm25_total := pmax(0, pm25_nrd) + pmax(0, pm25_onr)]

  # save annual xy table
  out_tot_rds <- file.path(out_dir, sprintf("TOTAL_%d_annual_xy.rds", yy))
  saveRDS(tot[, .(x,y,pm25_total,pm25_nrd,pm25_onr)], out_tot_rds)

  # ---- build coarse raster from table (anchored) ----
  r_tot <- annual_xy_to_raster(tot[, .(x,y,pm25_total)], ex_road, cmaq_res, p4s)
  names(r_tot) <- sprintf("AUTO_TOTAL_annual_%d", yy)

  out_coarse <- file.path(out_dir, sprintf("AUTO_TOTAL_annual_%d_anchored.tif", yy))
  writeRaster(r_tot, out_coarse, overwrite=TRUE)

  # ---- mean-preserving downscale to 1km ----
  tot_on_1km <- resample(r_tot, road_r, method="ngb")
  road_mean_on_1km <- resample(road_mean_12km, road_r, method="ngb")

  eps <- 1e-12
  den <- road_mean_on_1km
  den[is.na(den) | den < eps] <- NA

  factor_mp <- road_r / den
  factor_mp[is.na(factor_mp) | !is.finite(factor_mp)] <- 1

  down_1km <- tot_on_1km * factor_mp
  names(down_1km) <- sprintf("DOWN1KM_AUTO_TOTAL_%s_%d", tag, yy)

  out_down <- file.path(out_dir, sprintf("DOWN1KM_AUTO_TOTAL_%s_%d.tif", tag, yy))
  writeRaster(down_1km, out_down, overwrite=TRUE)

  # ---- QC summaries ----
  sum_tot <- tot[, .(
    n_cells = .N,
    min  = min(pm25_total, na.rm=TRUE),
    p50  = quantile(pm25_total, 0.50, na.rm=TRUE),
    p95  = quantile(pm25_total, 0.95, na.rm=TRUE),
    p99  = quantile(pm25_total, 0.99, na.rm=TRUE),
    max  = max(pm25_total, na.rm=TRUE),
    mean = mean(pm25_total, na.rm=TRUE),
    sd   = sd(pm25_total, na.rm=TRUE),
    share_gt200 = mean(pm25_total > 200, na.rm=TRUE)
  )]

  m0 <- cellStats(tot_on_1km, "mean")
  m1 <- cellStats(down_1km, "mean")
  mx <- cellStats(down_1km, "max")
  mn <- cellStats(down_1km, "min")

  # aggregate back to 12km and compare (mean-preservation check)
  down_back_12km <- aggregate(down_1km, fact=fact, fun=mean, na.rm=TRUE)
  down_back_12km <- resample(down_back_12km, r_tot, method="bilinear")
  diff_r <- down_back_12km - r_tot

  absdiff <- cellStats(abs(diff_r), "max")
  meandiff <- cellStats(diff_r, "mean")

  cat("Annual TOTAL (table) summary:\n",
      paste(capture.output(print(sum_tot)), collapse="\n"),
      "\n\n", file=qc_report, append=TRUE)

  cat(sprintf("Downscale QC:\n  mean_coarse=%.6f | mean_down1km=%.6f\n  min_down1km=%.6f | max_down1km=%.6f\n  max_absdiff_back12km=%.6f\n  mean_diff_back12km=%.6f\n\n",
              cellStats(r_tot,"mean"), m1, mn, mx, absdiff, meandiff),
      file=qc_report, append=TRUE)

  cat("Outputs:\n",
      "  ", out_coarse, "\n",
      "  ", out_down, "\n",
      "  ", out_tot_rds, "\n",
      "  ", qc_report, "\n\n", file=qc_report, append=TRUE)

  cat("Saved:\n  ", out_coarse, "\n  ", out_down, "\n")
  rm(nrd, onr, pm_nrd, pm_onr, tot, r_tot, tot_on_1km, down_1km, down_back_12km, diff_r); gc()
}

cat("\nDONE ✅\nOutputs written to:\n  ", out_dir, "\n")
