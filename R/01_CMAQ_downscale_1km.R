# ============================================================
# CMAQ annual mean (2011-2020) + optional downscale to roadiness
# Author: (for Sherry)
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(stringr)
  library(data.table)
})

# -----------------------------
# 0) USER SETTINGS
# -----------------------------
post_dir <- "/home/xshan2/HAQ_LAB/Sumaiya/cmaq/cmaq_output/POST"

# Choose which CMAQ post-processed file family to use
# Typical: hr2day_*.nc (daily mean surface) or COMBINE_ACONC_*.nc (hourly, huge)
# Strongly recommended: hr2day
prefer_pattern <- "hr2day"     # change if your filenames differ

# Target variable name inside netCDF (try this first)
var_prefer <- "PM25_TOT"       # you can change to PM25, PM25_TOT_ONR, etc.

# Output folders
out_dir_annual <- file.path(post_dir, "ANNUAL_PM25")
dir.create(out_dir_annual, showWarnings = FALSE, recursive = TRUE)

# Optional: downscale to roadiness (set to FALSE if you only want annual CMAQ)
do_downscale <- TRUE

# Your roadiness raster paths (choose 1km or 500m; you can run twice)
roadiness_path <- "/home/xshan2/your_path/roadiness_1km.tif"   # <-- CHANGE THIS
# roadiness_path <- "/home/xshan2/your_path/roadiness_500m.tif" # <-- OR THIS

out_dir_downscaled <- file.path(post_dir, "DOWNSCALED_PM25")
dir.create(out_dir_downscaled, showWarnings = FALSE, recursive = TRUE)

# Years to process
years <- 2011:2020

# -----------------------------
# 1) HELPERS
# -----------------------------

# Try to open a netCDF variable as a SpatRaster
open_cmaq_var <- function(nc_file, var_prefer = "PM25_TOT") {
  # Some netCDFs contain subdatasets; terra::sds() lists them
  s <- try(terra::sds(nc_file), silent = TRUE)

  if (!inherits(s, "try-error")) {
    subnames <- names(s)

    # 1) exact match
    if (var_prefer %in% subnames) return(rast(s[var_prefer]))

    # 2) heuristic match: contains PM25 and TOT or just PM25
    cand <- subnames[str_detect(subnames, regex("PM\\s*25|PM25", ignore_case = TRUE))]
    if (length(cand) > 0) {
      # prefer something like TOT
      cand2 <- cand[str_detect(cand, regex("TOT", ignore_case = TRUE))]
      pick <- if (length(cand2) > 0) cand2[1] else cand[1]
      message("  [var auto-pick] Using variable: ", pick)
      return(rast(s[pick]))
    }

    # 3) If none found, show variables
    stop("No PM2.5-like variable found. Available variables:\n  - ",
         paste(subnames, collapse = "\n  - "))
  }

  # Fallback: try direct open (works if file is simple)
  r <- try(rast(nc_file, subds = var_prefer), silent = TRUE)
  if (!inherits(r, "try-error")) return(r)

  stop("Could not open netCDF with terra. File: ", nc_file)
}

# Parse year-month from filename (robust-ish)
# Looks for YYYYMM or YYYY-MM in filename
parse_ym_from_filename <- function(f) {
  bn <- basename(f)
  # try YYYYMM
  m1 <- str_match(bn, "(20\\d{2})(0[1-9]|1[0-2])")
  if (!is.na(m1[1,1])) return(list(year = as.integer(m1[1,2]), month = as.integer(m1[1,3])))

  # try YYYY-MM
  m2 <- str_match(bn, "(20\\d{2})[-_](0[1-9]|1[0-2])")
  if (!is.na(m2[1,1])) return(list(year = as.integer(m2[1,2]), month = as.integer(m2[1,3])))

  # try YYYY_Mon (rare)
  stop("Cannot parse year/month from filename: ", bn)
}

days_in_month <- function(year, month) {
  # month 1..12
  start <- as.Date(sprintf("%04d-%02d-01", year, month))
  end <- if (month == 12) as.Date(sprintf("%04d-01-01", year + 1)) else as.Date(sprintf("%04d-%02d-01", year, month + 1))
  as.integer(end - start)
}

# Annual mean from monthly files (weighted by days)
compute_annual_from_monthlies <- function(files, year, var_prefer) {
  if (length(files) == 0) stop("No monthly files found for year ", year)

  sum_r <- NULL
  total_days <- 0L

  for (f in files) {
    ym <- parse_ym_from_filename(f)
    nd <- days_in_month(ym$year, ym$month)

    message(sprintf("  -> %s (days=%d)", basename(f), nd))
    r <- open_cmaq_var(f, var_prefer = var_prefer)

    # If hr2day: usually many layers (days). We need monthly mean across time layers
    # terra::app streams; no need to load entire cube
    r_month_mean <- app(r, fun = mean, na.rm = TRUE)

    if (is.null(sum_r)) {
      sum_r <- r_month_mean * nd
    } else {
      # Align if needed
      if (!compareGeom(sum_r, r_month_mean, stopOnError = FALSE)) {
        r_month_mean <- resample(r_month_mean, sum_r, method = "near")
      }
      sum_r <- sum_r + (r_month_mean * nd)
    }
    total_days <- total_days + nd
  }

  annual <- sum_r / total_days
  names(annual) <- paste0("PM25_annual_", year)
  annual
}

# Mean-preserving downscale to a finer roadiness grid (approx mass/mean preserved per CMAQ cell)
downscale_mean_preserving <- function(cmaq_annual, road_fine) {
  # Reproject roadiness to CMAQ CRS (or CMAQ to roadiness). We do everything in roadiness CRS to output fine grid.
  # Step 1: project CMAQ annual to road CRS but keep coarse cell values (nearest)
  cmaq_in_road_crs <- project(cmaq_annual, crs(road_fine), method = "near")

  # Step 2: assign each fine cell the coarse CMAQ value (nearest)
  cmaq_on_fine <- resample(cmaq_in_road_crs, road_fine, method = "near")

  # Step 3: compute coarse-cell mean roadiness (on the projected CMAQ grid)
  road_on_coarse <- resample(road_fine, cmaq_in_road_crs, method = "bilinear")  # average-ish proxy

  # Avoid division by zero / NA: set very small epsilon
  eps <- 1e-12
  road_on_coarse[is.na(road_on_coarse)] <- 0
  road_fine2 <- road_fine
  road_fine2[is.na(road_fine2)] <- 0

  # Step 4: bring coarse mean roadiness back to fine grid (nearest)
  road_mean_on_fine <- resample(road_on_coarse, road_fine2, method = "near")

  # Step 5: mean-preserving factor within each coarse cell:
  # factor = road_fine / mean_road_coarse
  # when mean_road_coarse == 0, set factor = 1 (no modulation)
  denom <- road_mean_on_fine
  denom[denom < eps] <- NA
  factor <- road_fine2 / denom
  factor[is.na(factor)] <- 1

  # Step 6: downscaled concentration
  down <- cmaq_on_fine * factor
  down
}

# -----------------------------
# 2) FIND CMAQ FILES
# -----------------------------
all_nc <- list.files(post_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)

# Keep only preferred family (e.g., hr2day)
cands <- all_nc[str_detect(tolower(basename(all_nc)), tolower(prefer_pattern))]

if (length(cands) == 0) {
  stop("No netCDF files found matching pattern '", prefer_pattern, "' under:\n  ", post_dir)
}

# Build table with parsed year/month
meta <- rbindlist(lapply(cands, function(f) {
  ym <- try(parse_ym_from_filename(f), silent = TRUE)
  if (inherits(ym, "try-error")) return(NULL)
  data.table(file = f, year = ym$year, month = ym$month)
}), fill = TRUE)

meta <- meta[year %in% years]
setorder(meta, year, month)

if (nrow(meta) == 0) stop("No files after filtering to years: ", paste(years, collapse = ","))

message("Found files by year:\n", capture.output(print(meta[, .N, by = year])))

# -----------------------------
# 3) COMPUTE ANNUAL CMAQ
# -----------------------------
annual_files <- list()

for (yy in years) {
  message("\n==============================")
  message("Computing annual mean for year: ", yy)
  message("==============================")

  files_y <- meta[year == yy][order(month)]$file
  if (length(files_y) == 0) {
    message("  [skip] No files for ", yy)
    next
  }

  r_annual <- compute_annual_from_monthlies(files_y, year = yy, var_prefer = var_prefer)

  out_tif <- file.path(out_dir_annual, sprintf("CMAQ_%s_annual_%d.tif", var_prefer, yy))
  writeRaster(r_annual, out_tif, overwrite = TRUE)
  message("  [saved] ", out_tif)

  annual_files[[as.character(yy)]] <- out_tif
}

# -----------------------------
# 4) OPTIONAL: DOWNSCALE TO ROADINESS (1km or 500m)
# -----------------------------
if (do_downscale) {
  message("\nLoading roadiness raster: ", roadiness_path)
  road <- rast(roadiness_path)

  for (yy in names(annual_files)) {
    message("\n------------------------------")
    message("Downscaling year: ", yy)
    message("------------------------------")

    cmaq_annual <- rast(annual_files[[yy]])

    down <- downscale_mean_preserving(cmaq_annual, road)
    names(down) <- paste0("PM25_downscaled_", yy)

    out_tif2 <- file.path(out_dir_downscaled, sprintf("CMAQ_%s_downscaled_%s.tif", var_prefer, yy))
    writeRaster(down, out_tif2, overwrite = TRUE)
    message("  [saved] ", out_tif2)
  }
}

message("\nDONE.")
