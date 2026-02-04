# ============================================================
# CMAQ (hr2day_w_organics) Annual (2011–2020) + 1km Downscaling
# Variables: PM25_TOT_ONR, PM25_TOT_NRD
# Dependencies: raster, ncdf4, data.table
# Output: GeoTIFFs in /scratch/xshan2/R_Code/disperseR/Auto
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(raster)
  library(ncdf4)
  library(data.table)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
years <- 2011:2020

# CMAQ hr2day_w_organics directory (copied to your scratch)
post_dir <- "/scratch/xshan2/cmaq_POST"

# Output directory
out_dir <- "/scratch/xshan2/R_Code/disperseR/Auto"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Projection (match your roadiness)
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000 +units=m +no_defs"

# CMAQ grid resolution (12 km)
cmaq_res <- 12000

# Roadiness CSV (1km)
road_csv <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/roadiness_1km_hw_loc_sherry.csv"

# Choose which normalized roadiness weight to use
road_weight_col <- "leng.distm2_hw_norm"   # highway
# road_weight_col <- "leng.distm2_loc_norm"  # local

# Tag used in output filenames
tag <- "hw"  # "hw" or "loc"

# CMAQ variables to process
vars <- c("PM25_TOT_ONR", "PM25_TOT_NRD")

# Save 12km annual rasters too?
save_coarse_annual <- TRUE

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

# ============================================================
# 0) READ ROADINESS -> 1km RASTER TEMPLATE
# ============================================================
cat("Reading roadiness CSV:\n  ", road_csv, "\n")

road_dt <- fread(road_csv)
cn <- names(road_dt)

xcol <- if ("x" %in% cn) "x" else if ("X" %in% cn) "X" else stop("No x column in roadiness CSV")
ycol <- if ("y" %in% cn) "y" else if ("Y" %in% cn) "Y" else stop("No y column in roadiness CSV")
stopifnot(road_weight_col %in% cn)

# Keep only what we need (reduces memory)
road_dt <- road_dt[, .(x = get(xcol), y = get(ycol), w = get(road_weight_col))]
road_dt[!is.finite(w) | is.na(w), w := 0]

road_r <- rasterFromXYZ(as.data.frame(road_dt), crs = CRS(p4s))
names(road_r) <- paste0("road_", tag)

cat("\n[Roadiness raster template]\n")
print(road_r)

# Aggregation factor (12km / 1km = 12)
fact <- round(cmaq_res / res(road_r)[1])
if (!is.finite(fact) || fact < 1) stop("Bad aggregation factor. Check roadiness resolution.")
cat("\nAggregation factor (coarse/fine): ", fact, "\n")

# Precompute 12km mean roadiness on the roadiness-anchored grid
road_mean_12km <- aggregate(road_r, fact = fact, fun = mean, na.rm = TRUE)

# Anchoring: we will build all CMAQ rasters using roadiness xmin/ymin as the origin
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

# De-duplicate: if a year-month has multiple files, keep the first by filename
meta <- meta[, .SD[1], by = .(year, month)]
setorder(meta, year, month)

cat("\nFiles used per year (after de-dup by year-month):\n")
print(meta[, .N, by = year])

# ============================================================
# 2) MAIN LOOP: YEARS x VARS
# ============================================================
for (yy in years) {

  files_y <- meta[year == yy]
  if (nrow(files_y) == 0) {
    cat("\n[SKIP] No files found for year:", yy, "\n")
    next
  }

  if (nrow(files_y) < 12) {
    warning("Year ", yy, " has only ", nrow(files_y), " months (expected 12). Continuing anyway.")
  }

  cat("\n========================================\n")
  cat("YEAR:", yy, "| months:", paste(sprintf("%02d", files_y$month), collapse = ","), "\n")
  cat("========================================\n")

  for (v in vars) {

    cat("\n  ---- Processing variable:", v, "----\n")

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

      # Monthly mean over time dimension
      pm_mean <- apply(pm, c(1, 2), mean, na.rm = TRUE)

      if (is.null(sum_mat)) {
        sum_mat <- pm_mean * nd
        ncol_cmaq <- dim(pm_mean)[1]
        nrow_cmaq <- dim(pm_mean)[2]
      } else {
        sum_mat <- sum_mat + pm_mean * nd
      }

      tot_days <- tot_days + nd
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
    # Step 1: assign each 1km cell the parent 12km value (nearest neighbor)
    cmaq_on_1km <- resample(r_cmaq, road_r, method = "ngb")

    # Step 2: compute mean roadiness within each 12km cell, then map back to 1km
    road_mean_on_1km <- resample(road_mean_12km, road_r, method = "ngb")

    # Step 3: factor = road / mean(road within 12km); handle zeros safely
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
  }
}

cat("\nDONE ✅\nOutputs written to:\n  ", out_dir, "\n")
