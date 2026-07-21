#!/usr/bin/env Rscript

# ==============================================================================
# Figure 4-5: AADT modification of wind-adjusted mobile-source PM2.5
# Virginia, 2011-2020
#
# Confirmed source definitions used throughout this script:
#   ONR = on-road + non-road gasoline sources
#   NRD = on-road + non-road diesel sources
#
# This script intentionally starts from the validated annual output of the
# Figure 4-4 wind model. It does not reimplement or replace that wind model.
#
# For 1-km cell i in parent CMAQ cell g and source s:
#
#   A*_i       = (AADT_i + a0) / mean_g(AADT + a0)
#   Q_i        = W_wind_i * (A*_i ^ alpha)
#   W_final_i  = Q_i / mean_g(Q)
#   PM_final_i = PM12_g * W_final_i
#
# The parent-cell arithmetic mean of PM_final therefore equals PM12_g.
#
# Default AADT mapping:
#   ONR (gasoline): aadt_total
#   NRD (diesel):   aadt_total
#
# This common total-AADT modifier matches the earlier analysis and is complete
# for 2011-2020. A source-specific sensitivity option (truck AADT for NRD) is
# available below, but it is rejected automatically when any annual truck field
# is absent or entirely zero.
#
# Scientific interpretation:
# AADT is a road-traffic activity proxy. Because the CMAQ source tags also
# include non-road gasoline/diesel emissions, AADT should be described as a
# spatial modifier of the combined mobile-source fields, not as a direct
# inventory of non-road activity.
#
# Main outputs:
#   Figure_4_5_wind_AADT_VA_2011_2020.pdf   (vector master)
#   Figure_4_5_wind_AADT_VA_2011_2020.tiff  (600 dpi, when supported)
#   Figure_4_5_wind_AADT_VA_2011_2020.png   (600 dpi, when supported)
#   Figure_4_5_annual_alpha_*.fst
#   Figure_4_5_decade_alpha_*.fst
#   Figure_4_5_year_QC.csv
#   Figure_4_5_parent_QC.csv
#   Figure_4_5_AADT_assignment_QC.csv
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fst)
  library(ggplot2)
  library(USAboundaries)
  library(scales)
  library(grid)
})

run_figure_4_5 <- function() {

SCRIPT_VERSION <- "Figure 4-5 wind + real annual HPMS AADT, strict 2011-2020 v3"

if (interactive()) {
  stop(
    "Interactive R session detected. Do not paste this script at the `>` ",
    "prompt. Exit R with q(save = \"no\") and run the complete file with: ",
    "Rscript --vanilla Figure_4_5_Wind_AADT_2011_2020_Nature.R"
  )
}

run_started_at <- Sys.time()

cat(
  "\n============================================================\n",
  SCRIPT_VERSION,
  "\n============================================================\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------------------------

years <- 2011:2020

p4s <- paste(
  "+proj=lcc",
  "+lat_1=33",
  "+lat_2=45",
  "+lat_0=40",
  "+lon_0=-97",
  "+a=6370000",
  "+b=6370000"
)

base_dir <- "/scratch/xshan2/R_Code/Automobiles"

figure44_dir <- file.path(
  base_dir,
  "FIGURES",
  "NATIVE_CMAQ_REBUILT",
  "FIGURE_4_4_WIND_DIRECTION_SPEED_FINAL_STRICT_2011_2020"
)

wind_plot_data_file <- file.path(
  figure44_dir,
  "Figure_4_4_plot_data.rds"
)

hpms_dir <- file.path(
  base_dir,
  "HPMS_AADT"
)

out_dir <- file.path(
  base_dir,
  "FIGURES",
  "NATIVE_CMAQ_REBUILT",
  "FIGURE_4_5_WIND_AADT_REAL_HPMS_2011_2020"
)

cache_dir <- file.path(
  out_dir,
  "AADT_GRID_CACHE"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Rebuild only the new caches in cache_dir. This never reads the old
# AADT_nearest_2020.fst produced by a carry-forward workflow.
force_rebuild_aadt_cache <- FALSE

# Sensitivity values. alpha = 1 is the main analysis and figure.
alpha_values <- c(0.5, 1, 2)
main_alpha <- 1

# AADT mapping used in the primary analysis.
#
# "total_for_both" (default and recommended for the current 2011-2020 files):
#   total AADT is a common road-activity modifier for both combined mobile tags.
#
# "total_gasoline_truck_diesel" (optional sensitivity):
#   total AADT modifies ONR and truck AADT modifies NRD. The script stops if the
#   truck field is missing/all-zero in even one year; it never silently falls
#   back to wind-only or mixes definitions across years.
aadt_mapping_mode <- "total_for_both"

valid_aadt_mapping_modes <- c(
  "total_for_both",
  "total_gasoline_truck_diesel"
)

if (!(aadt_mapping_mode %in% valid_aadt_mapping_modes)) {
  stop(
    "Unknown aadt_mapping_mode: ",
    aadt_mapping_mode,
    ". Valid choices: ",
    paste(valid_aadt_mapping_modes, collapse = ", ")
  )
}

aadt_mapping_description <- if (
  aadt_mapping_mode == "total_for_both"
) {
  paste(
    "Annual total HPMS AADT was used as a common road-traffic activity",
    "modifier for both gasoline and diesel mobile-source fields."
  )
} else {
  paste(
    "Annual total HPMS AADT was used for gasoline and annual truck HPMS",
    "AADT was used for diesel."
  )
}

# One vehicle day-1 is negligible relative to observed AADT, but prevents
# zero-AADT cells from receiving an exactly zero mathematical weight.
aadt_offset <- 1

# Strict numerical tolerance for parent-cell mean conservation.
conservation_tolerance <- 1e-8

# Main-map display limits use a common quantile cap across both sources.
# Values above the cap are squished, not removed. The cap is written to CSV.
absolute_upper_quantile <- 0.995
difference_upper_quantile <- 0.995

# Figure dimensions and output resolution.
figure_width_in <- 12.2
figure_height_in <- 7.0
raster_dpi <- 600

# Annual 2 x 5 maps are useful as supplementary figures but are not recommended
# as the main six-panel figure.
make_annual_supplement <- TRUE

pdf_file <- file.path(
  out_dir,
  "Figure_4_5_wind_AADT_VA_2011_2020.pdf"
)

tiff_file <- file.path(
  out_dir,
  "Figure_4_5_wind_AADT_VA_2011_2020.tiff"
)

png_file <- file.path(
  out_dir,
  "Figure_4_5_wind_AADT_VA_2011_2020.png"
)

year_qc_file <- file.path(out_dir, "Figure_4_5_year_QC.csv")
parent_qc_file <- file.path(out_dir, "Figure_4_5_parent_QC.csv")
aadt_qc_file <- file.path(out_dir, "Figure_4_5_AADT_assignment_QC.csv")
scale_file <- file.path(out_dir, "Figure_4_5_display_scales.csv")
caption_file <- file.path(out_dir, "Figure_4_5_caption.txt")
method_file <- file.path(out_dir, "Figure_4_5_method.txt")
completion_file <- file.path(out_dir, "Figure_4_5_RUN_COMPLETE.txt")

# A completion marker is valid only for the current successful run.
if (file.exists(completion_file)) {
  unlink(completion_file)
}

alpha_tag <- function(alpha) {
  gsub(
    "\\.",
    "p",
    format(alpha, scientific = FALSE, trim = TRUE)
  )
}

annual_output_file <- function(alpha) {
  file.path(
    out_dir,
    sprintf(
      "Figure_4_5_annual_alpha_%s.fst",
      alpha_tag(alpha)
    )
  )
}

decade_output_file <- function(alpha) {
  file.path(
    out_dir,
    sprintf(
      "Figure_4_5_decade_alpha_%s.fst",
      alpha_tag(alpha)
    )
  )
}

annual_total_output_file <- function(alpha) {
  file.path(
    out_dir,
    sprintf(
      "Figure_4_5_annual_TOTAL_alpha_%s.fst",
      alpha_tag(alpha)
    )
  )
}

decade_total_output_file <- function(alpha) {
  file.path(
    out_dir,
    sprintf(
      "Figure_4_5_decade_TOTAL_alpha_%s.fst",
      alpha_tag(alpha)
    )
  )
}

# ------------------------------------------------------------------------------
# 2. Strict preflight
# ------------------------------------------------------------------------------

if (!file.exists(wind_plot_data_file)) {
  stop(
    "Missing validated Figure 4-4 data:\n",
    wind_plot_data_file
  )
}

hpms_files <- setNames(
  file.path(
    hpms_dir,
    sprintf(
      "Virginia_%d_HPMS_AADT_sf.rds",
      years
    )
  ),
  years
)

missing_hpms <- hpms_files[!file.exists(hpms_files)]

if (length(missing_hpms) > 0) {
  stop(
    "Missing annual HPMS files:\n",
    paste(missing_hpms, collapse = "\n")
  )
}

wind_object <- readRDS(wind_plot_data_file)

if (!is.list(wind_object) || !("annual_data" %in% names(wind_object))) {
  stop(
    "Figure_4_4_plot_data.rds does not contain annual_data."
  )
}

annual_wind <- as.data.table(wind_object$annual_data)

required_wind_columns <- c(
  "x",
  "y",
  "parent_id",
  "year",
  "source_code",
  "pm12",
  "pm_wind",
  "wind_weight",
  "raw_wind_weight"
)

missing_wind_columns <- setdiff(
  required_wind_columns,
  names(annual_wind)
)

if (length(missing_wind_columns) > 0) {
  stop(
    "annual_data is missing required columns: ",
    paste(missing_wind_columns, collapse = ", ")
  )
}

annual_wind <- annual_wind[
  year %in% years &
    source_code %in% c("ONR", "NRD")
]

annual_wind[, `:=`(
  x = as.numeric(x),
  y = as.numeric(y),
  parent_id = as.integer(parent_id),
  year = as.integer(year),
  pm12 = as.numeric(pm12),
  pm_wind = as.numeric(pm_wind),
  wind_weight = as.numeric(wind_weight),
  raw_wind_weight = as.numeric(raw_wind_weight)
)]

if (
  anyDuplicated(
    annual_wind[, .(x, y, year, source_code)]
  ) > 0
) {
  stop(
    "Duplicate x/y/year/source records were found in Figure 4-4 annual_data."
  )
}

expected_source_years <- CJ(
  source_code = c("NRD", "ONR"),
  year = years
)

actual_source_years <- unique(
  annual_wind[, .(source_code, year)]
)

setkey(expected_source_years, source_code, year)
setkey(actual_source_years, source_code, year)

missing_source_years <- fsetdiff(
  expected_source_years,
  actual_source_years
)

if (nrow(missing_source_years) > 0) {
  stop(
    "Figure 4-4 annual_data is incomplete:\n",
    paste(capture.output(print(missing_source_years)), collapse = "\n")
  )
}

annual_wind[, source_name := fifelse(
  source_code == "ONR",
  "Gasoline mobile sources",
  "Diesel mobile sources"
)]

nonfinite_wind <- annual_wind[
  !is.finite(x) |
    !is.finite(y) |
    !is.finite(parent_id) |
    !is.finite(pm12) |
    !is.finite(pm_wind) |
    !is.finite(wind_weight) |
    !is.finite(raw_wind_weight)
]

if (nrow(nonfinite_wind) > 0) {
  stop(
    "Non-finite values were found in validated wind annual_data."
  )
}

wind_parent_check <- annual_wind[, .(
  pm12 = first(pm12),
  mean_pm_wind = mean(pm_wind),
  mean_wind_weight = mean(wind_weight)
), by = .(year, source_code, parent_id)]

wind_parent_check[, wind_error := mean_pm_wind - pm12]

maximum_existing_wind_error <- max(
  abs(wind_parent_check$wind_error),
  na.rm = TRUE
)

if (maximum_existing_wind_error > conservation_tolerance) {
  stop(
    "The Figure 4-4 wind input fails parent conservation. Maximum error = ",
    signif(maximum_existing_wind_error, 8)
  )
}

# A single unique 1-km grid is required for every source and year.
grid_dt <- unique(
  annual_wind[, .(x, y, parent_id)]
)

if (anyDuplicated(grid_dt[, .(x, y)]) > 0) {
  stop(
    "The same 1-km coordinate maps to more than one parent_id."
  )
}

expected_rows <- nrow(grid_dt) * length(years) * 2L

if (nrow(annual_wind) != expected_rows) {
  stop(
    "Unexpected annual_data size. Expected ",
    expected_rows,
    " rows but found ",
    nrow(annual_wind),
    "."
  )
}

cat(
  "Preflight passed: ",
  nrow(grid_dt),
  " 1-km cells; both sources; all years 2011-2020.\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 3. Build or read annual 1-km AADT assignments
# ------------------------------------------------------------------------------

read_hpms_year <- function(year) {

  file_path <- unname(
    hpms_files[as.character(year)]
  )

  hpms <- readRDS(file_path)

  if (!inherits(hpms, "sf")) {
    stop(
      "HPMS annual RDS is not an sf object: ",
      file_path
    )
  }

  required_fields <- c(
    "aadt_total",
    "aadt_truck"
  )

  missing_fields <- setdiff(
    required_fields,
    names(hpms)
  )

  if (length(missing_fields) > 0) {
    stop(
      "HPMS file is missing standardized fields: ",
      file_path,
      "\nMissing: ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (is.na(st_crs(hpms))) {
    stop(
      "HPMS file has no CRS: ",
      file_path
    )
  }

  hpms <- st_zm(
    hpms,
    drop = TRUE,
    what = "ZM"
  )

  hpms <- hpms[
    !st_is_empty(hpms),
    c("aadt_total", "aadt_truck")
  ]

  hpms$aadt_total <- suppressWarnings(
    as.numeric(hpms$aadt_total)
  )

  hpms$aadt_truck <- suppressWarnings(
    as.numeric(hpms$aadt_truck)
  )

  hpms$aadt_total[
    !is.finite(hpms$aadt_total) |
      hpms$aadt_total < 0
  ] <- 0

  hpms$aadt_truck[
    !is.finite(hpms$aadt_truck) |
      hpms$aadt_truck < 0
  ] <- 0

  hpms <- st_transform(
    hpms,
    st_crs(p4s)
  )

  hpms
}

build_aadt_grid_year <- function(year) {

  cache_file <- file.path(
    cache_dir,
    sprintf(
      "VA_1km_real_HPMS_AADT_%d.fst",
      year
    )
  )

  if (
    file.exists(cache_file) &&
      !force_rebuild_aadt_cache
  ) {

    cached <- read_fst(
      cache_file,
      as.data.table = TRUE
    )

    required_cache_columns <- c(
      "x",
      "y",
      "year",
      "aadt_total",
      "aadt_truck"
    )

    if (
      all(required_cache_columns %in% names(cached)) &&
        nrow(cached) == nrow(grid_dt) &&
        identical(
          cached[order(x, y), .(x, y)],
          grid_dt[order(x, y), .(x, y)]
        )
    ) {
      cat(
        "Using validated AADT cache for ",
        year,
        ": ",
        cache_file,
        "\n",
        sep = ""
      )

      return(cached)
    }

    warning(
      "Ignoring invalid AADT cache and rebuilding: ",
      cache_file
    )
  }

  cat(
    "Reading real HPMS AADT for ",
    year,
    " ...\n",
    sep = ""
  )

  hpms <- read_hpms_year(year)

  if (nrow(hpms) == 0) {
    stop(
      "HPMS has zero usable features for ",
      year,
      "."
    )
  }

  points_sf <- st_as_sf(
    grid_dt,
    coords = c("x", "y"),
    crs = st_crs(p4s),
    remove = FALSE
  )

  # st_nearest_feature uses a spatial index. One call per year avoids rebuilding
  # that index for every source or small chunk.
  nearest_index <- st_nearest_feature(
    points_sf,
    hpms
  )

  if (
    length(nearest_index) != nrow(grid_dt) ||
      anyNA(nearest_index)
  ) {
    stop(
      "Nearest-feature assignment failed for HPMS ",
      year,
      "."
    )
  }

  out <- copy(grid_dt[, .(x, y)])
  out[, year := as.integer(year)]
  out[, aadt_total := as.numeric(hpms$aadt_total[nearest_index])]
  out[, aadt_truck := as.numeric(hpms$aadt_truck[nearest_index])]

  out[
    !is.finite(aadt_total) |
      aadt_total < 0,
    aadt_total := 0
  ]

  out[
    !is.finite(aadt_truck) |
      aadt_truck < 0,
    aadt_truck := 0
  ]

  write_fst(
    out,
    cache_file
  )

  cat(
    "WROTE AADT cache: ",
    cache_file,
    " | rows=",
    nrow(out),
    "\n",
    sep = ""
  )

  out
}

aadt_grid_list <- vector(
  "list",
  length(years)
)

names(aadt_grid_list) <- as.character(years)

for (year in years) {
  aadt_grid_list[[as.character(year)]] <- build_aadt_grid_year(year)
  gc()
}

aadt_grid <- rbindlist(
  aadt_grid_list,
  use.names = TRUE
)

setkey(aadt_grid, year, x, y)
setkey(annual_wind, year, x, y)

annual_base <- aadt_grid[annual_wind]

if (
  anyNA(annual_base$aadt_total) ||
    anyNA(annual_base$aadt_truck)
) {
  stop(
    "Some Figure 4-4 cells did not receive annual AADT values."
  )
}

if (aadt_mapping_mode == "total_for_both") {
  annual_base[, aadt_raw := aadt_total]
} else {
  annual_base[, aadt_raw := fifelse(
    source_code == "ONR",
    aadt_total,
    aadt_truck
  )]
}

aadt_assignment_qc <- annual_base[, .(
  n_1km_cells = .N,
  n_parent_cells = uniqueN(parent_id),
  aadt_min = min(aadt_raw),
  aadt_mean = mean(aadt_raw),
  aadt_median = median(aadt_raw),
  aadt_p95 = as.numeric(quantile(aadt_raw, 0.95)),
  aadt_p99 = as.numeric(quantile(aadt_raw, 0.99)),
  aadt_max = max(aadt_raw),
  percent_zero = 100 * mean(aadt_raw == 0)
), by = .(year, source_code, source_name)]

fwrite(
  aadt_assignment_qc,
  aadt_qc_file
)

cat(
  "\nAADT assignment summary:\n",
  sep = ""
)

print(aadt_assignment_qc)

if (aadt_mapping_mode == "total_gasoline_truck_diesel") {
  invalid_truck_years <- aadt_assignment_qc[
    source_code == "NRD" &
      (!is.finite(aadt_max) | aadt_max <= 0),
    year
  ]

  if (length(invalid_truck_years) > 0) {
    stop(
      "Truck AADT is absent or entirely zero for: ",
      paste(invalid_truck_years, collapse = ", "),
      ". Do not mix wind-only NRD years with truck-AADT NRD years. ",
      "Repair the historical truck fields or use ",
      "aadt_mapping_mode = \"total_for_both\"."
    )
  }
}

# ------------------------------------------------------------------------------
# 4. Combine wind and AADT; normalize within every parent CMAQ cell
# ------------------------------------------------------------------------------

calculate_wind_aadt <- function(base_dt, aadt_alpha) {

  out <- copy(base_dt)

  out[, aadt_parent_mean := mean(
    aadt_raw + aadt_offset
  ), by = .(year, source_code, parent_id)]

  out[, aadt_relative :=
        (aadt_raw + aadt_offset) /
          aadt_parent_mean]

  out[, raw_wind_aadt_weight :=
        wind_weight *
          (aadt_relative ^ aadt_alpha)]

  out[, parent_raw_weight_mean := mean(
    raw_wind_aadt_weight
  ), by = .(year, source_code, parent_id)]

  out[, wind_aadt_weight := fifelse(
    is.finite(parent_raw_weight_mean) &
      parent_raw_weight_mean > 0,
    raw_wind_aadt_weight /
      parent_raw_weight_mean,
    wind_weight
  )]

  if (
    any(!is.finite(out$wind_aadt_weight)) ||
      any(out$wind_aadt_weight < 0)
  ) {
    stop(
      "Invalid final weights for alpha=",
      aadt_alpha,
      "."
    )
  }

  out[, pm_wind_aadt := pm12 * wind_aadt_weight]
  out[, delta_aadt := pm_wind_aadt - pm_wind]

  # Figure 4-4 annual_data may itself contain a column named `alpha`. Using a
  # different function argument and the `..` prefix prevents data.table from
  # reading that inherited column instead of the current AADT sensitivity.
  out[, alpha := as.numeric(..aadt_alpha)]

  out
}

year_qc_list <- list()
parent_qc_list <- list()
main_annual <- NULL
main_decade <- NULL

for (aadt_alpha in alpha_values) {

  cat(
    "\nProcessing alpha=",
    aadt_alpha,
    " ...\n",
    sep = ""
  )

  result <- calculate_wind_aadt(
    annual_base,
    aadt_alpha
  )

  parent_qc <- result[, .(
    pm12 = first(pm12),
    mean_pm_wind = mean(pm_wind),
    mean_pm_wind_aadt = mean(pm_wind_aadt),
    mean_wind_weight = mean(wind_weight),
    mean_wind_aadt_weight = mean(wind_aadt_weight),
    n_1km_cells = .N
  ), by = .(alpha, year, source_code, source_name, parent_id)]

  parent_qc[, wind_error := mean_pm_wind - pm12]
  parent_qc[, wind_aadt_error := mean_pm_wind_aadt - pm12]
  parent_qc[, wind_aadt_weight_error := mean_wind_aadt_weight - 1]

  maximum_error <- max(
    abs(parent_qc$wind_aadt_error),
    na.rm = TRUE
  )

  if (maximum_error > conservation_tolerance) {
    stop(
      "Parent conservation failed for alpha=",
      aadt_alpha,
      ". Maximum absolute error=",
      signif(maximum_error, 8)
    )
  }

  year_qc <- result[, .(
    n_1km_cells = .N,
    n_parent_cells = uniqueN(parent_id),
    mean_pm_wind = mean(pm_wind),
    mean_pm_wind_aadt = mean(pm_wind_aadt),
    mean_delta = mean(delta_aadt),
    mean_absolute_delta = mean(abs(delta_aadt)),
    median_absolute_delta = median(abs(delta_aadt)),
    p95_absolute_delta = as.numeric(
      quantile(abs(delta_aadt), 0.95)
    ),
    delta_p05 = as.numeric(quantile(delta_aadt, 0.05)),
    delta_p95 = as.numeric(quantile(delta_aadt, 0.95)),
    percent_absolute_change_above_0_01 =
      100 * mean(abs(delta_aadt) >= 0.01),
    percent_absolute_change_above_0_05 =
      100 * mean(abs(delta_aadt) >= 0.05)
  ), by = .(alpha, year, source_code, source_name)]

  year_conservation_qc <- parent_qc[, .(
    maximum_absolute_conservation_error = max(
      abs(wind_aadt_error),
      na.rm = TRUE
    )
  ), by = .(alpha, year, source_code, source_name)]

  year_qc <- merge(
    year_qc,
    year_conservation_qc,
    by = c(
      "alpha",
      "year",
      "source_code",
      "source_name"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  decade <- result[, .(
    pm12 = mean(pm12),
    pm_wind = mean(pm_wind),
    pm_wind_aadt = mean(pm_wind_aadt),
    delta_aadt = mean(delta_aadt),
    mean_aadt = mean(aadt_raw),
    mean_aadt_relative = mean(aadt_relative),
    mean_wind_weight = mean(wind_weight),
    mean_wind_aadt_weight = mean(wind_aadt_weight),
    n_years = uniqueN(year)
  ), by = .(
    x,
    y,
    parent_id,
    source_code,
    source_name,
    alpha
  )]

  if (any(decade$n_years != length(years))) {
    stop(
      "Decade output is missing one or more years for alpha=",
      aadt_alpha,
      ". n_years distribution: ",
      paste(
        names(table(decade$n_years)),
        as.integer(table(decade$n_years)),
        sep = ":",
        collapse = ", "
      )
    )
  }

  annual_to_write <- result[, .(
    x,
    y,
    parent_id,
    year,
    source_code,
    source_name,
    alpha,
    pm12,
    pm_wind,
    aadt_total,
    aadt_truck,
    aadt_raw,
    aadt_relative,
    wind_weight,
    wind_aadt_weight,
    pm_wind_aadt,
    delta_aadt
  )]

  gasoline_total_part <- result[
    source_code == "ONR",
    .(
      x,
      y,
      parent_id,
      year,
      alpha,
      pm12_gasoline = pm12,
      pm_wind_gasoline = pm_wind,
      pm_wind_aadt_gasoline = pm_wind_aadt,
      delta_aadt_gasoline = delta_aadt
    )
  ]

  diesel_total_part <- result[
    source_code == "NRD",
    .(
      x,
      y,
      parent_id,
      year,
      alpha,
      pm12_diesel = pm12,
      pm_wind_diesel = pm_wind,
      pm_wind_aadt_diesel = pm_wind_aadt,
      delta_aadt_diesel = delta_aadt
    )
  ]

  total_annual <- merge(
    gasoline_total_part,
    diesel_total_part,
    by = c(
      "x",
      "y",
      "parent_id",
      "year",
      "alpha"
    ),
    all = FALSE,
    sort = FALSE
  )

  if (
    nrow(total_annual) !=
      nrow(grid_dt) * length(years)
  ) {
    stop(
      "Gasoline and diesel annual grids did not align for TOTAL, alpha=",
      aadt_alpha,
      "."
    )
  }

  total_annual[, `:=`(
    source_code = "TOTAL",
    source_name = "Total gasoline + diesel mobile sources",
    pm12 = pm12_gasoline + pm12_diesel,
    pm_wind = pm_wind_gasoline + pm_wind_diesel,
    pm_wind_aadt =
      pm_wind_aadt_gasoline +
        pm_wind_aadt_diesel,
    delta_aadt =
      delta_aadt_gasoline +
        delta_aadt_diesel
  )]

  total_decade <- total_annual[, .(
    pm12 = mean(pm12),
    pm_wind = mean(pm_wind),
    pm_wind_aadt = mean(pm_wind_aadt),
    delta_aadt = mean(delta_aadt),
    n_years = uniqueN(year)
  ), by = .(
    x,
    y,
    parent_id,
    source_code,
    source_name,
    alpha
  )]

  total_annual_to_write <- total_annual[, .(
    x,
    y,
    parent_id,
    year,
    source_code,
    source_name,
    alpha,
    pm12,
    pm_wind,
    pm_wind_aadt,
    delta_aadt
  )]

  write_fst(
    annual_to_write,
    annual_output_file(aadt_alpha)
  )

  write_fst(
    decade,
    decade_output_file(aadt_alpha)
  )

  write_fst(
    total_annual_to_write,
    annual_total_output_file(aadt_alpha)
  )

  write_fst(
    total_decade,
    decade_total_output_file(aadt_alpha)
  )

  cat(
    "WROTE: ",
    annual_output_file(aadt_alpha),
    "\nWROTE: ",
    decade_output_file(aadt_alpha),
    "\nWROTE: ",
    annual_total_output_file(aadt_alpha),
    "\nWROTE: ",
    decade_total_output_file(aadt_alpha),
    "\n",
    sep = ""
  )

  year_qc_list[[alpha_tag(aadt_alpha)]] <- year_qc
  parent_qc_list[[alpha_tag(aadt_alpha)]] <- parent_qc

  if (isTRUE(all.equal(aadt_alpha, main_alpha))) {
    main_annual <- annual_to_write
    main_decade <- decade
  }

  rm(
    result,
    parent_qc,
    year_qc,
    year_conservation_qc,
    decade,
    annual_to_write,
    gasoline_total_part,
    diesel_total_part,
    total_annual,
    total_decade,
    total_annual_to_write
  )

  gc()
}

if (is.null(main_annual) || is.null(main_decade)) {
  stop(
    "main_alpha was not included in alpha_values."
  )
}

year_qc_all <- rbindlist(
  year_qc_list,
  use.names = TRUE
)

parent_qc_all <- rbindlist(
  parent_qc_list,
  use.names = TRUE
)

fwrite(year_qc_all, year_qc_file)
fwrite(parent_qc_all, parent_qc_file)

# ------------------------------------------------------------------------------
# 5. Publication figure: 2011-2020 mean, wind versus wind + AADT
# ------------------------------------------------------------------------------

states_ll <- USAboundaries::us_states(
  resolution = "low"
)

virginia_lcc <- states_ll[
  states_ll$state_abbr == "VA",
]

virginia_lcc <- st_transform(
  st_as_sf(virginia_lcc),
  st_crs(p4s)
)

va_bbox <- st_bbox(virginia_lcc)
va_xlim <- c(va_bbox[["xmin"]], va_bbox[["xmax"]])
va_ylim <- c(va_bbox[["ymin"]], va_bbox[["ymax"]])

absolute_values <- c(
  main_decade$pm_wind,
  main_decade$pm_wind_aadt
)

absolute_limit <- as.numeric(
  quantile(
    absolute_values,
    absolute_upper_quantile,
    na.rm = TRUE
  )
)

if (!is.finite(absolute_limit) || absolute_limit <= 0) {
  absolute_limit <- max(absolute_values, na.rm = TRUE)
}

difference_limit <- as.numeric(
  quantile(
    abs(main_decade$delta_aadt),
    difference_upper_quantile,
    na.rm = TRUE
  )
)

if (!is.finite(difference_limit) || difference_limit <= 0) {
  difference_limit <- max(
    abs(main_decade$delta_aadt),
    na.rm = TRUE
  )
}

absolute_breaks <- pretty(
  c(0, absolute_limit),
  n = 5
)

absolute_breaks <- absolute_breaks[
  absolute_breaks >= 0 &
    absolute_breaks <= absolute_limit
]

difference_breaks <- pretty(
  c(-difference_limit, difference_limit),
  n = 5
)

difference_breaks <- difference_breaks[
  difference_breaks >= -difference_limit &
    difference_breaks <= difference_limit
]

fwrite(
  data.table(
    scale = c("absolute", "difference"),
    lower = c(0, -difference_limit),
    upper = c(absolute_limit, difference_limit),
    quantile_cap = c(
      absolute_upper_quantile,
      difference_upper_quantile
    ),
    alpha = main_alpha
  ),
  scale_file
)

absolute_colours <- c(
  "#FFFFE5",
  "#FEE391",
  "#FEC44F",
  "#FE9929",
  "#D95F0E",
  "#993404"
)

difference_colours <- c(
  "#2166AC",
  "#67A9CF",
  "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7",
  "#EF8A62",
  "#B2182B"
)

common_map_layers <- function() {
  list(
    geom_sf(
      data = virginia_lcc,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey20",
      linewidth = 0.28
    ),
    coord_sf(
      crs = st_crs(p4s),
      xlim = va_xlim,
      ylim = va_ylim,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ),
    theme_void(
      base_size = 9,
      base_family = "sans"
    ),
    theme(
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.border = element_rect(
        fill = NA,
        color = "grey45",
        linewidth = 0.25
      ),
      plot.title = element_text(
        size = 9.5,
        face = "bold",
        hjust = 0,
        margin = margin(b = 3)
      ),
      plot.margin = margin(
        t = 3,
        r = 4,
        b = 2,
        l = 4,
        unit = "pt"
      )
    )
  )
}

make_absolute_panel <- function(
    data,
    value_column,
    panel_title,
    show_legend = FALSE) {

  plot_dt <- copy(data)
  plot_dt[, value_to_plot := as.numeric(get(value_column))]

  ggplot() +
    geom_raster(
      data = plot_dt,
      aes(x = x, y = y, fill = value_to_plot),
      interpolate = FALSE
    ) +
    common_map_layers() +
    scale_fill_gradientn(
      colours = absolute_colours,
      limits = c(0, absolute_limit),
      breaks = absolute_breaks,
      labels = label_number(accuracy = 0.1),
      oob = squish,
      na.value = "white",
      name = expression(PM[2.5]~(mu*g~m^{-3})),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = unit(76, "mm"),
        barheight = unit(3.5, "mm"),
        frame.colour = "grey35",
        frame.linewidth = 0.25,
        ticks.colour = "grey35"
      )
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = NULL
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 8.8),
      legend.text = element_text(size = 8.0)
    )
}

make_difference_panel <- function(
    data,
    panel_title,
    show_legend = FALSE) {

  ggplot() +
    geom_raster(
      data = data,
      aes(x = x, y = y, fill = delta_aadt),
      interpolate = FALSE
    ) +
    common_map_layers() +
    scale_fill_gradientn(
      colours = difference_colours,
      values = rescale(
        c(
          -difference_limit,
          -difference_limit * 0.66,
          -difference_limit * 0.33,
          0,
          difference_limit * 0.33,
          difference_limit * 0.66,
          difference_limit
        )
      ),
      limits = c(-difference_limit, difference_limit),
      breaks = difference_breaks,
      labels = label_number(
        accuracy = if (difference_limit < 0.1) 0.01 else 0.1
      ),
      oob = squish,
      na.value = "white",
      name = expression(Delta*PM[2.5]~(mu*g~m^{-3})),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = unit(46, "mm"),
        barheight = unit(3.5, "mm"),
        frame.colour = "grey35",
        frame.linewidth = 0.25,
        ticks.colour = "grey35"
      )
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = NULL
    ) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.title = element_text(size = 8.8),
      legend.text = element_text(size = 8.0)
    )
}

gasoline_data <- main_decade[source_code == "ONR"]
diesel_data <- main_decade[source_code == "NRD"]

p_a <- make_absolute_panel(
  gasoline_data,
  "pm_wind",
  "(a)  Gasoline: wind",
  FALSE
)

p_b <- make_absolute_panel(
  gasoline_data,
  "pm_wind_aadt",
  "(b)  Gasoline: wind + AADT",
  FALSE
)

p_c <- make_difference_panel(
  gasoline_data,
  "(c)  Gasoline: AADT effect",
  FALSE
)

p_d <- make_absolute_panel(
  diesel_data,
  "pm_wind",
  "(d)  Diesel: wind",
  FALSE
)

p_e <- make_absolute_panel(
  diesel_data,
  "pm_wind_aadt",
  "(e)  Diesel: wind + AADT",
  FALSE
)

p_f <- make_difference_panel(
  diesel_data,
  "(f)  Diesel: AADT effect",
  FALSE
)

get_legend <- function(plot_object) {

  grob <- ggplotGrob(plot_object)
  guide_ids <- grep("^guide-box", grob$layout$name)

  for (i in guide_ids) {
    candidate <- grob$grobs[[i]]
    if (!inherits(candidate, "zeroGrob")) {
      return(candidate)
    }
  }

  stop("No non-empty legend was found.")
}

absolute_legend <- get_legend(
  make_absolute_panel(
    gasoline_data,
    "pm_wind_aadt",
    "",
    TRUE
  ) +
    theme(
      panel.border = element_blank(),
      plot.title = element_blank()
    )
)

difference_legend <- get_legend(
  make_difference_panel(
    gasoline_data,
    "",
    TRUE
  ) +
    theme(
      panel.border = element_blank(),
      plot.title = element_blank()
    )
)

panel_grobs <- lapply(
  list(p_a, p_b, p_c, p_d, p_e, p_f),
  ggplotGrob
)

draw_main_figure <- function() {

  grid.newpage()

  layout <- grid.layout(
    nrow = 3,
    ncol = 3,
    widths = unit(c(1, 1, 1), "null"),
    heights = unit.c(
      unit(1, "null"),
      unit(1, "null"),
      unit(0.68, "in")
    )
  )

  pushViewport(viewport(layout = layout))

  draw_at <- function(grob, row, column) {
    pushViewport(
      viewport(
        layout.pos.row = row,
        layout.pos.col = column
      )
    )
    grid.draw(grob)
    popViewport()
  }

  draw_at(panel_grobs[[1]], 1, 1)
  draw_at(panel_grobs[[2]], 1, 2)
  draw_at(panel_grobs[[3]], 1, 3)
  draw_at(panel_grobs[[4]], 2, 1)
  draw_at(panel_grobs[[5]], 2, 2)
  draw_at(panel_grobs[[6]], 2, 3)

  pushViewport(
    viewport(
      layout.pos.row = 3,
      layout.pos.col = 1:2
    )
  )
  grid.draw(absolute_legend)
  popViewport()

  pushViewport(
    viewport(
      layout.pos.row = 3,
      layout.pos.col = 3
    )
  )
  grid.draw(difference_legend)
  popViewport()

  popViewport()
}

pdf_partial_file <- paste0(pdf_file, ".partial")

pdf(
  pdf_partial_file,
  width = figure_width_in,
  height = figure_height_in,
  useDingbats = FALSE,
  family = "sans"
)
draw_main_figure()
dev.off()

if (
  !file.exists(pdf_partial_file) ||
    file.info(pdf_partial_file)$size < 10000
) {
  stop("The vector figure was not rendered correctly: ", pdf_partial_file)
}

if (!file.copy(pdf_partial_file, pdf_file, overwrite = TRUE)) {
  stop("Could not install the completed vector figure: ", pdf_file)
}

unlink(pdf_partial_file)

cat(
  "WROTE vector figure: ",
  pdf_file,
  "\n",
  sep = ""
)

save_raster_figure <- function(file_path, file_type) {

  ok <- FALSE

  if (requireNamespace("ragg", quietly = TRUE)) {
    tryCatch({
      if (file_type == "png") {
        ragg::agg_png(
          file_path,
          width = figure_width_in,
          height = figure_height_in,
          units = "in",
          res = raster_dpi,
          background = "white"
        )
      } else {
        ragg::agg_tiff(
          file_path,
          width = figure_width_in,
          height = figure_height_in,
          units = "in",
          res = raster_dpi,
          compression = "lzw",
          background = "white"
        )
      }

      draw_main_figure()
      dev.off()
      ok <- file.exists(file_path) && file.info(file_path)$size > 0
    }, error = function(e) {
      warning(
        "ragg could not create ",
        file_type,
        ": ",
        conditionMessage(e)
      )
      while (dev.cur() > 1) dev.off()
    })
  }

  if (!ok && capabilities("cairo")) {
    tryCatch({
      if (file_type == "png") {
        png(
          file_path,
          width = figure_width_in,
          height = figure_height_in,
          units = "in",
          res = raster_dpi,
          type = "cairo-png",
          bg = "white"
        )
      } else {
        tiff(
          file_path,
          width = figure_width_in,
          height = figure_height_in,
          units = "in",
          res = raster_dpi,
          compression = "lzw",
          type = "cairo",
          bg = "white"
        )
      }

      draw_main_figure()
      dev.off()
      ok <- file.exists(file_path) && file.info(file_path)$size > 0
    }, error = function(e) {
      warning(
        "Cairo could not create ",
        file_type,
        ": ",
        conditionMessage(e)
      )
      while (dev.cur() > 1) dev.off()
    })
  }

  ok
}

png_created <- save_raster_figure(
  png_file,
  "png"
)

tiff_created <- save_raster_figure(
  tiff_file,
  "tiff"
)

# Headless HPC installations sometimes lack ragg/Cairo devices. In that case,
# convert the already validated vector PDF with Poppler when pdftocairo exists.
convert_pdf_with_poppler <- function(file_path, file_type) {

  converter <- Sys.which("pdftocairo")

  if (!nzchar(converter)) {
    return(FALSE)
  }

  output_prefix <- sub(
    "\\.(png|tif|tiff)$",
    "",
    file_path,
    ignore.case = TRUE
  )

  format_flag <- if (file_type == "png") "-png" else "-tiff"

  conversion_log <- suppressWarnings(
    system2(
      converter,
      args = c(
        format_flag,
        "-singlefile",
        "-r",
        as.character(raster_dpi),
        pdf_file,
        output_prefix
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  )

  exit_status <- attr(conversion_log, "status")
  if (is.null(exit_status)) exit_status <- 0L

  created_file <- paste0(
    output_prefix,
    if (file_type == "png") ".png" else ".tif"
  )

  if (
    exit_status != 0L ||
      !file.exists(created_file) ||
      file.info(created_file)$size <= 0
  ) {
    warning(
      "pdftocairo could not create ",
      file_type,
      ": ",
      paste(conversion_log, collapse = " | ")
    )
    return(FALSE)
  }

  if (!identical(created_file, file_path)) {
    if (!file.copy(created_file, file_path, overwrite = TRUE)) {
      return(FALSE)
    }
    unlink(created_file)
  }

  file.exists(file_path) && file.info(file_path)$size > 0
}

if (!png_created) {
  png_created <- convert_pdf_with_poppler(png_file, "png")
}

if (!tiff_created) {
  tiff_created <- convert_pdf_with_poppler(tiff_file, "tiff")
}

# ------------------------------------------------------------------------------
# 6. Annual supplementary maps
# ------------------------------------------------------------------------------

if (make_annual_supplement) {

  annual_plot <- as.data.table(copy(main_annual))
  annual_plot[, year_f := factor(year, levels = years)]

  for (source in c("ONR", "NRD")) {

    source_label <- if (source == "ONR") {
      "Gasoline mobile sources"
    } else {
      "Diesel mobile sources"
    }

    source_dt <- annual_plot[source_code == source]

    p_annual <- ggplot() +
      geom_raster(
        data = source_dt,
        aes(x = x, y = y, fill = pm_wind_aadt),
        interpolate = FALSE
      ) +
      geom_sf(
        data = virginia_lcc,
        inherit.aes = FALSE,
        fill = NA,
        color = "grey25",
        linewidth = 0.18
      ) +
      coord_sf(
        crs = st_crs(p4s),
        xlim = va_xlim,
        ylim = va_ylim,
        expand = FALSE,
        datum = NA
      ) +
      facet_wrap(
        ~year_f,
        ncol = 5
      ) +
      scale_fill_gradientn(
        colours = absolute_colours,
        limits = c(0, absolute_limit),
        breaks = absolute_breaks,
        labels = label_number(accuracy = 0.1),
        oob = squish,
        na.value = "white",
        name = expression(PM[2.5]~(mu*g~m^{-3})),
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(65, "mm"),
          barheight = unit(3.5, "mm")
        )
      ) +
      labs(
        title = paste0(
          source_label,
          ": wind + AADT (alpha = ",
          main_alpha,
          ")"
        ),
        x = NULL,
        y = NULL
      ) +
      theme_void(
        base_size = 9,
        base_family = "sans"
      ) +
      theme(
        strip.text = element_text(
          size = 8.5,
          face = "bold"
        ),
        plot.title = element_text(
          size = 10,
          face = "bold",
          hjust = 0
        ),
        legend.position = "bottom",
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        panel.spacing = unit(2, "mm"),
        plot.margin = margin(4, 4, 4, 4)
      )

    annual_pdf <- file.path(
      out_dir,
      sprintf(
        "Figure_S_annual_%s_wind_AADT_2011_2020.pdf",
        source
      )
    )

    pdf(
      annual_pdf,
      width = 12.2,
      height = 5.7,
      useDingbats = FALSE,
      family = "sans"
    )
    print(p_annual)
    dev.off()

    cat(
      "WROTE annual supplement: ",
      annual_pdf,
      "\n",
      sep = ""
    )
  }
}

# ------------------------------------------------------------------------------
# 7. Caption, method record and completion marker
# ------------------------------------------------------------------------------

caption_text <- paste0(
  "Figure 4-5 | Influence of traffic activity on wind-adjusted, source-specific ",
  "PM2.5 in Virginia. Maps show 2011-2020 mean PM2.5 attributed to gasoline ",
  "mobile sources (a-c) and diesel mobile sources (d-f) after downscaling with ",
  "roadiness and annual wind fields (a,d), and after additionally incorporating ",
  "annual HPMS AADT (b,e). Difference maps (c,f) show wind + AADT minus wind. ",
  aadt_mapping_description,
  " Weights were ",
  "normalized within each parent CMAQ cell, preserving its annual mean ",
  "source-specific concentration. Display limits correspond to the ",
  absolute_upper_quantile * 100,
  "th percentile for absolute concentrations and the ",
  difference_upper_quantile * 100,
  "th percentile of absolute differences; values beyond these limits are ",
  "shown at the endpoint colours."
)

writeLines(
  caption_text,
  caption_file
)

method_text <- c(
  paste0("Script version: ", SCRIPT_VERSION),
  "Source definitions: ONR = on-road + non-road gasoline; NRD = on-road + non-road diesel.",
  "Wind input: validated annual_data from Figure 4-4.",
  "AADT input: real annual Virginia HPMS sf RDS files for every year 2011-2020.",
  paste0("AADT mapping mode: ", aadt_mapping_mode, "."),
  aadt_mapping_description,
  paste0("Alpha values: ", paste(alpha_values, collapse = ", "), "; main alpha: ", main_alpha, "."),
  paste0("AADT offset: ", aadt_offset, " vehicle day-1."),
  "Spatial assignment: value of the nearest annual HPMS road feature at each 1-km cell centre.",
  "Normalization: combined wind-AADT weights have parent-cell arithmetic mean one.",
  "Interpretation: AADT is a road-traffic proxy applied to combined on-road + non-road source tags.",
  paste0("Maximum existing wind conservation error: ", signif(maximum_existing_wind_error, 8)),
  paste0("Vector PDF: ", pdf_file),
  paste0("PNG created: ", png_created),
  paste0("TIFF created: ", tiff_created)
)

writeLines(
  method_text,
  method_file
)

finalize_run <- function() {

  alpha_output_files <- unlist(
    lapply(alpha_values, function(aadt_alpha) {
      c(
        annual_output_file(aadt_alpha),
        decade_output_file(aadt_alpha),
        annual_total_output_file(aadt_alpha),
        decade_total_output_file(aadt_alpha)
      )
    }),
    use.names = FALSE
  )

  supplementary_files <- character()

  if (make_annual_supplement) {
    supplementary_files <- file.path(
      out_dir,
      sprintf(
        "Figure_S_annual_%s_wind_AADT_2011_2020.pdf",
        c("ONR", "NRD")
      )
    )
  }

  required_current_outputs <- c(
    pdf_file,
    year_qc_file,
    parent_qc_file,
    aadt_qc_file,
    scale_file,
    caption_file,
    method_file,
    alpha_output_files,
    supplementary_files
  )

  output_exists <- file.exists(required_current_outputs)
  output_info <- file.info(required_current_outputs)
  output_is_current <- output_exists &
    is.finite(output_info$size) &
    output_info$size > 0 &
    output_info$mtime >= (run_started_at - 2)

  if (any(!output_is_current)) {
    stop(
      "Run is incomplete; no completion marker was written. Missing, empty, ",
      "or stale outputs:\n",
      paste(
        required_current_outputs[!output_is_current],
        collapse = "\n"
      )
    )
  }

  writeLines(
    c(
      paste0("Completed: ", Sys.time()),
      paste0("Script: ", SCRIPT_VERSION),
      "Years: 2011-2020",
      "Sources: ONR gasoline; NRD diesel",
      paste0("AADT mapping: ", aadt_mapping_mode),
      paste0("Main alpha: ", main_alpha),
      paste0("PDF: ", pdf_file),
      paste0("PNG created: ", png_created),
      paste0("TIFF created: ", tiff_created)
    ),
    completion_file
  )

  cat(
    "\n============================================================\n",
    "FIGURE 4-5 COMPLETE\n",
    "Output directory: ",
    out_dir,
    "\nPDF: ",
    pdf_file,
    "\nPNG created: ",
    png_created,
    "\nTIFF created: ",
    tiff_created,
    "\n============================================================\n",
    sep = ""
  )
}

finalize_run()

}

run_figure_4_5()
