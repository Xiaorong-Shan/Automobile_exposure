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
read_year_file <- function(dir_path, pattern_year) {
  # e.g., pattern_year="2019-01_2019-12"
  f <- file.path(dir_path, sprintf("PM25_TOT_%s_cmaq_%s.rds",
                                   basename(dir_path), pattern_year))
  # basename(dir_path) is "NRD"/"ONR" -> but file naming uses NRD/ONR
  # safer: build explicitly below instead of this helper
  f
}

annual_xy_from_daily <- function(dt, valcol) {
  setDT(dt)
  dt[, .(
    pm25_mean = mean(get(valcol), na.rm = TRUE),
    pm25_p99  = quantile(get(valcol), 0.99, na.rm = TRUE),
    n_days    = sum(is.finite(get(valcol)))
  ), by = .(x, y)]
}

qc_sample_daily <- function(dt, valcol) {
  setDT(dt)
  set.seed(1)
  s <- dt[sample.int(.N, min(sample_n, .N)), get(valcol)]
  s <- s[is.finite(s)]
  list(
    q = quantile(s, c(0,.01,.05,.5,.95,.99,1), na.rm=TRUE),
    share_neg = mean(s < 0, na.rm=TRUE),
    share_gt50 = mean(s > 50, na.rm=TRUE),
    share_gt200 = mean(s > 200, na.rm=TRUE)
  )
}

build_coarse_raster_from_xy <- function(tot_xy, ex_road, cmaq_res, p4s) {
  # tot_xy: data.table with x,y,pm25_total on regular grid
  setDT(tot_xy)

  ux <- sort(unique(tot_xy$x))
  uy <- sort(unique(tot_xy$y))
  ncol_cmaq <- length(ux)
  nrow_cmaq <- length(uy)

  # long -> matrix (rows=y, cols=x)
  mdt <- dcast(tot_xy, y ~ x, value.var = "pm25_total")
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

downscale_mean_preserving <- function(r_coarse, road_r, road_mean_12km, fact, tag, year, out_dir) {
  # 12km -> 1km nearest
  coarse_on_1km <- resample(r_coarse, road_r, method = "ngb")

  # 12km mean roadiness -> 1km nearest
  road_mean_on_1km <- resample(road_mean_12km, road_r, method = "ngb")

  # factor = road / mean_road within coarse cell
  eps <- 1e-12
  den <- road_mean_on_1km
  den[is.na(den) | den < eps] <- NA

  factor_mp <- road_r / den
  factor_mp[is.na(factor_mp) | !is.finite(factor_mp)] <- 1

  down_1km <- coarse_on_1km * factor_mp
  names(down_1km) <- sprintf("AUTO_TOTAL_down1km_%s_%d", tag, year)

  out_down <- file.path(out_dir, sprintf("DOWN1KM_AUTO_TOTAL_%s_%d.tif", tag, year))
  writeRaster(down_1km, out_down, overwrite = TRUE)

  # QC: aggregate back
  down_back_12km <- aggregate(down_1km, fact = fact, fun = mean, na.rm = TRUE)
  down_back_12km <- resample(down_back_12km, r_coarse, method = "bilinear")

  diff_r <- down_back_12km - r_coarse

  list(
    out_down = out_down,
    mean_coarse = cellStats(r_coarse, "mean"),
    mean_down1km = cellStats(down_1km, "mean"),
    max_down1km = cellStats(down_1km, "max"),
    max_absdiff = cellStats(abs(diff_r), "max"),
    mean_diff = cellStats(diff_r, "mean")
  )
}

# ============================================================
# 0) ROADINESS TEMPLATE + 12km MEAN ROADINESS
# ============================================================
cat("Reading roadiness CSV:\n  ", road_csv, "\n")
road_dt <- fread(road_csv)
cn <- names(road_dt)

xcol <- if ("x" %in% cn) "x" else if ("X" %in% cn) "X" else stop("No x column in roadiness CSV")
ycol <- if ("y" %in% cn) "y" else if ("Y" %in% cn) "Y" else stop("No y column in roadiness CSV")
stopifnot(road_weight_col %in% cn)

road_dt <- road_dt[, .(x = as.numeric(get(xcol)),
                       y = as.numeric(get(ycol)),
                       w = as.numeric(get(road_weight_col)))]
road_dt[!is.finite(w) | is.na(w) | w < 0, w := 0]

road_r <- rasterFromXYZ(as.data.frame(road_dt), crs = CRS(p4s))
names(road_r) <- paste0("road_", tag)

fact <- round(cmaq_res / res(road_r)[1])
if (!is.finite(fact) || fact < 1) stop("Bad aggregation factor. Check roadiness resolution.")
road_mean_12km <- aggregate(road_r, fact = fact, fun = mean, na.rm = TRUE)
ex_road <- extent(road_r)

cat("\n[Roadiness]\n"); print(road_r)
cat("Aggregation factor (12km/1km): ", fact, "\n")

# ============================================================
# 1) LOOP YEARS
# ============================================================
for (yy in years) {

  cat("\n========================================\n")
  cat("YEAR:", yy, "\n")
  cat("========================================\n")

  # input file names (match your copied rds naming)
  pat <- sprintf("%d-01_%d-12", yy, yy)

  f_nrd <- file.path(nrd_dir, sprintf("PM25_TOT_NRD_cmaq_%s.rds", pat))
  f_onr <- file.path(onr_dir, sprintf("PM25_TOT_ONR_cmaq_%s.rds", pat))

  if (!file.exists(f_nrd)) {
    cat("[SKIP] Missing NRD:", f_nrd, "\n")
    next
  }
  if (!file.exists(f_onr)) {
    cat("[SKIP] Missing ONR:", f_onr, "\n")
    next
  }

  # ---- read daily long tables ----
  nrd <- readRDS(f_nrd); setDT(nrd)
  onr <- readRDS(f_onr); setDT(onr)

  # ---- quick daily QC (optional, lightweight) ----
  qc_nrd <- qc_sample_daily(nrd, "PM25_TOT_NRD")
  qc_onr <- qc_sample_daily(onr, "PM25_TOT_ONR")

  # ---- annual mean by (x,y) ----
  pm_nrd <- annual_xy_from_daily(nrd, "PM25_TOT_NRD")
  pm_onr <- annual_xy_from_daily(onr, "PM25_TOT_ONR")
  rm(nrd, onr); gc()

  # physical constraint
  pm_nrd[pm25_mean < 0 | !is.finite(pm25_mean), pm25_mean := 0]
  pm_onr[pm25_mean < 0 | !is.finite(pm25_mean), pm25_mean := 0]

  # ---- TOTAL ----
  setkey(pm_nrd, x, y)
  setkey(pm_onr, x, y)

  tot <- pm_nrd[pm_onr, on=.(x,y)]
  setnames(tot, c("pm25_mean", "i.pm25_mean"), c("pm25_nrd", "pm25_onr"))
  tot[, pm25_total := pm25_nrd + pm25_onr]
  tot[pm25_total < 0 | !is.finite(pm25_total), pm25_total := 0]

  # save annual table
  out_tot_dt <- file.path(out_dir, sprintf("TOTAL_%d_annual_xy.rds", yy))
  saveRDS(tot[, .(x, y, pm25_total, pm25_nrd, pm25_onr, n_days, i.n_days)], out_tot_dt)

  # ---- coarse raster (12km) ----
  r_tot <- build_coarse_raster_from_xy(tot[, .(x,y,pm25_total)], ex_road, cmaq_res, p4s)
  names(r_tot) <- sprintf("AUTO_TOTAL_annual_%d", yy)

  out_coarse <- file.path(out_dir, sprintf("AUTO_TOTAL_annual_%d_anchored.tif", yy))
  writeRaster(r_tot, out_coarse, overwrite = TRUE)

  # ---- downscale 1km ----
  qc_ds <- downscale_mean_preserving(r_tot, road_r, road_mean_12km, fact, tag, yy, out_dir)

  # ---- write QC txt ----
  qc_txt <- file.path(out_dir, sprintf("QC_AUTO_TOTAL_%s_%d.txt", tag, yy))
  writeLines(c(
    paste0("QC Auto TOTAL annual + 1km downscale (", tag, ") year ", yy),
    paste0("NRD rds: ", f_nrd),
    paste0("ONR rds: ", f_onr),
    paste0("Annual XY rds: ", out_tot_dt),
    paste0("Coarse tif: ", out_coarse),
    paste0("Down1km tif: ", qc_ds$out_down),
    "",
    "Daily sample QC (NRD):",
    paste0("  share_neg=", qc_nrd$share_neg, " | share_gt50=", qc_nrd$share_gt50, " | share_gt200=", qc_nrd$share_gt200),
    paste(capture.output(print(qc_nrd$q)), collapse="\n"),
    "",
    "Daily sample QC (ONR):",
    paste0("  share_neg=", qc_onr$share_neg, " | share_gt50=", qc_onr$share_gt50, " | share_gt200=", qc_onr$share_gt200),
    paste(capture.output(print(qc_onr$q)), collapse="\n"),
    "",
    "Downscale QC:",
    sprintf("  mean_coarse=%.6f | mean_down1km=%.6f", qc_ds$mean_coarse, qc_ds$mean_down1km),
    sprintf("  max_down1km=%.6f", qc_ds$max_down1km),
    sprintf("  max_absdiff_back12km=%.6f", qc_ds$max_absdiff),
    sprintf("  mean_diff_back12km=%.6f", qc_ds$mean_diff)
  ), qc_txt)

  cat("Saved:\n  ", out_coarse, "\n  ", qc_ds$out_down, "\n  ", out_tot_dt, "\n  ", qc_txt, "\n")
}

cat("\nDONE ✅\nOutputs in:\n  ", out_dir, "\n")
