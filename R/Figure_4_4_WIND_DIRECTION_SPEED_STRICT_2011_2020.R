#!/usr/bin/env Rscript

# ==============================================================================
# Figure 4-4
# Influence of prevailing wind direction and wind speed on roadiness-based
# gasoline- and diesel-attributable PM2.5 in Virginia, averaged over 2011-2020
#
# Confirmed source definitions:
#   ONR = gasoline vehicles
#   NRD = diesel vehicles
#
# Wind treatment (faithful to the earlier model formulation):
#   - Wind_E and Wind_N are used as signed vector components.
#   - Annual component vectors define prevailing transport direction.
#   - Eight neighboring 1-km roadiness values are weighted by angular alignment
#     with the upwind source direction: 1 / (absolute angular difference + 0.05).
#   - A relative wind-speed factor is applied:
#       ((domain mean wind speed + eps) / (local wind speed + U0))^beta_wind
#     with U0 = 0.1 m s-1 and beta_wind = 0.5.
#   - The final raw wind modifier is:
#       upwind roadiness x wind-speed factor.
#   - Weights are normalized within each parent CMAQ cell so that the arithmetic
#     mean of the 1-km concentrations equals the parent 12-km concentration.
#
# Figure panels:
#   (a) Gasoline roadiness-only
#   (b) Gasoline roadiness + wind
#   (c) Wind-adjusted minus roadiness-only, gasoline
#   (d) Diesel roadiness-only
#   (e) Diesel roadiness + wind
#   (f) Wind-adjusted minus roadiness-only, diesel
#
# Absolute concentration panels use 0-1.0 ug m-3.
# Difference panels use a common symmetric scale centered at zero.
#
# INPUTS:
#   Annual source-specific CMAQ means:
#     /scratch/xshan2/R_Code/Automobiles/FIGURES/NATIVE_CMAQ_REBUILT/
#       ANNUAL_ONR_mean_YYYY.rds
#       ANNUAL_NRD_mean_YYYY.rds
#
#   Figure 4-3 model-ready roadiness outputs:
#     .../FIGURE_4_3_ROADINESS/
#       VA_1km_gasoline_roadiness_2011_2020_mean.fst
#       VA_1km_diesel_roadiness_2011_2020_mean.fst
#
#   Annual NLDAS component fields:
#     /scratch/xshan2/NLDAS_wind/VA/
#       NLDAS_VA_YYYY_Wind_EN_mean.tif
#
# OUTPUTS:
#   Figure_4_4_wind_VA_2011_2020.pdf
#   Figure_4_4_wind_VA_2011_2020.png
#   Figure_4_4_QC_summary.csv
#   Figure_4_4_year_QC.csv
#   Figure_4_4_parent_QC.csv
#   Figure_4_4_plot_data.rds
#   VA_1km_gasoline_wind_2011_2020_mean.fst
#   VA_1km_diesel_wind_2011_2020_mean.fst
#
# HOPPER:
#   The script writes a vector PDF first and then tries multiple independent
#   headless PNG methods. It reports success only after all required outputs
#   actually exist.
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fst)
  library(raster)
  library(ggplot2)
  library(USAboundaries)
  library(scales)
  library(viridis)
})

SCRIPT_VERSION <- "Figure 4-4 direction+speed strict-v3 (complete 2011-2020 only)"

cat(
  "\n============================================================\n",
  SCRIPT_VERSION,
  "\n============================================================\n",
  "This script will not process any year unless all 2011-2020 inputs pass preflight.\n",
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

cmaq_dir <- paste0(
  "/scratch/xshan2/R_Code/Automobiles/FIGURES/",
  "NATIVE_CMAQ_REBUILT"
)

figure43_dir <- file.path(
  cmaq_dir,
  "FIGURE_4_3_ROADINESS"
)

gasoline_template_file <- file.path(
  figure43_dir,
  "VA_1km_gasoline_roadiness_2011_2020_mean.fst"
)

diesel_template_file <- file.path(
  figure43_dir,
  "VA_1km_diesel_roadiness_2011_2020_mean.fst"
)

gasoline_decade_file <- file.path(
  cmaq_dir,
  "DECADE_ONR_mean_2011_2020.rds"
)

diesel_decade_file <- file.path(
  cmaq_dir,
  "DECADE_NRD_mean_2011_2020.rds"
)

wind_search_root <- "/scratch/xshan2/NLDAS_wind"
preferred_wind_dir <- file.path(
  wind_search_root,
  "VA"
)

out_dir <- file.path(
  cmaq_dir,
  "FIGURE_4_4_WIND_DIRECTION_SPEED_FINAL_STRICT_2011_2020"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf_file <- file.path(
  out_dir,
  "Figure_4_4_wind_VA_2011_2020.pdf"
)

png_file <- file.path(
  out_dir,
  "Figure_4_4_wind_VA_2011_2020.png"
)

summary_qc_file <- file.path(
  out_dir,
  "Figure_4_4_QC_summary.csv"
)

year_qc_file <- file.path(
  out_dir,
  "Figure_4_4_year_QC.csv"
)

parent_qc_file <- file.path(
  out_dir,
  "Figure_4_4_parent_QC.csv"
)

plot_data_file <- file.path(
  out_dir,
  "Figure_4_4_plot_data.rds"
)

gasoline_output_file <- file.path(
  out_dir,
  "VA_1km_gasoline_wind_2011_2020_mean.fst"
)

diesel_output_file <- file.path(
  out_dir,
  "VA_1km_diesel_wind_2011_2020_mean.fst"
)

preflight_inventory_file <- file.path(
  out_dir,
  "Figure_4_4_preflight_inventory.csv"
)

completion_file <- file.path(
  out_dir,
  "Figure_4_4_RUN_COMPLETE.txt"
)

# A PNG is useful for Word/preview, but the vector PDF is the scientific master.
# Set TRUE only if this Hopper R environment has a working headless PNG renderer.
require_png <- FALSE

# ------------------------------------------------------------------------------
# 1A. Strict preflight: resolve and validate every 2011-2020 input
# ------------------------------------------------------------------------------

static_inputs <- data.table(
  input_type = c(
    "Figure 4-3 gasoline template",
    "Figure 4-3 diesel template",
    "Gasoline decade CMAQ",
    "Diesel decade CMAQ"
  ),
  year = NA_integer_,
  expected_file = c(
    gasoline_template_file,
    diesel_template_file,
    gasoline_decade_file,
    diesel_decade_file
  )
)

annual_cmaq_inputs <- rbindlist(
  list(
    data.table(
      input_type = "Annual gasoline CMAQ",
      year = years,
      expected_file = file.path(
        cmaq_dir,
        sprintf(
          "ANNUAL_ONR_mean_%d.rds",
          years
        )
      )
    ),
    data.table(
      input_type = "Annual diesel CMAQ",
      year = years,
      expected_file = file.path(
        cmaq_dir,
        sprintf(
          "ANNUAL_NRD_mean_%d.rds",
          years
        )
      )
    )
  ),
  use.names = TRUE
)

all_wind_tifs <- if (
  dir.exists(
    wind_search_root
  )
) {
  list.files(
    path = wind_search_root,
    pattern = "\\.tif$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
} else {
  character()
}

resolve_wind_file <- function(
    year) {

  expected_basename <- sprintf(
    "NLDAS_VA_%d_Wind_EN_mean.tif",
    year
  )

  preferred_file <- file.path(
    preferred_wind_dir,
    expected_basename
  )

  if (
    file.exists(
      preferred_file
    )
  ) {
    return(
      list(
        path = preferred_file,
        resolution_status = "preferred exact path",
        candidate_count = 1L
      )
    )
  }

  exact_hits <- all_wind_tifs[
    basename(
      all_wind_tifs
    ) ==
      expected_basename
  ]

  if (
    length(
      exact_hits
    ) ==
      1
  ) {
    return(
      list(
        path = exact_hits,
        resolution_status = "exact basename found recursively",
        candidate_count = 1L
      )
    )
  }

  if (
    length(
      exact_hits
    ) >
      1
  ) {
    return(
      list(
        path = paste(
          exact_hits,
          collapse = " | "
        ),
        resolution_status = "ambiguous exact basename",
        candidate_count = length(
          exact_hits
        )
      )
    )
  }

  relaxed_pattern <- paste0(
    "(?i)VA.*",
    year,
    ".*Wind[_-]?EN.*mean.*\\.tif$"
  )

  relaxed_hits <- all_wind_tifs[
    grepl(
      relaxed_pattern,
      basename(
        all_wind_tifs
      ),
      perl = TRUE
    )
  ]

  if (
    length(
      relaxed_hits
    ) ==
      1
  ) {
    return(
      list(
        path = relaxed_hits,
        resolution_status = "relaxed filename match found recursively",
        candidate_count = 1L
      )
    )
  }

  if (
    length(
      relaxed_hits
    ) >
      1
  ) {
    return(
      list(
        path = paste(
          relaxed_hits,
          collapse = " | "
        ),
        resolution_status = "ambiguous relaxed filename match",
        candidate_count = length(
          relaxed_hits
        )
      )
    )
  }

  list(
    path = NA_character_,
    resolution_status = "not found",
    candidate_count = 0L
  )
}

wind_resolution <- lapply(
  years,
  resolve_wind_file
)

wind_inventory <- rbindlist(
  lapply(
    seq_along(
      years
    ),
    function(
        i) {

      data.table(
        input_type = "Annual NLDAS wind",
        year = years[i],
        expected_file = file.path(
          preferred_wind_dir,
          sprintf(
            "NLDAS_VA_%d_Wind_EN_mean.tif",
            years[i]
          )
        ),
        resolved_file = wind_resolution[[i]]$path,
        resolution_status = wind_resolution[[i]]$resolution_status,
        candidate_count = wind_resolution[[i]]$candidate_count
      )
    }
  ),
  use.names = TRUE
)

static_inputs[, resolved_file := expected_file]
static_inputs[, resolution_status := fifelse(
  file.exists(
    expected_file
  ),
  "exact path",
  "not found"
)]
static_inputs[, candidate_count := as.integer(
  file.exists(
    expected_file
  )
)]

annual_cmaq_inputs[, resolved_file := expected_file]
annual_cmaq_inputs[, resolution_status := fifelse(
  file.exists(
    expected_file
  ),
  "exact path",
  "not found"
)]
annual_cmaq_inputs[, candidate_count := as.integer(
  file.exists(
    expected_file
  )
)]

preflight_inventory <- rbindlist(
  list(
    static_inputs,
    annual_cmaq_inputs,
    wind_inventory
  ),
  use.names = TRUE,
  fill = TRUE
)

preflight_inventory[, exists := !is.na(
  resolved_file
) &
  candidate_count ==
  1L &
  file.exists(
    resolved_file
  )]

preflight_inventory[, size_mb := fifelse(
  exists,
  round(
    file.info(
      resolved_file
    )$size /
      1024^2,
    3
  ),
  NA_real_
)]

# Validate each resolved wind file before any model year is processed.
preflight_inventory[, `:=`(
  wind_layers = NA_integer_,
  wind_layer_names = NA_character_,
  wind_valid = NA
)]

for (
  i in which(
    preflight_inventory$input_type ==
      "Annual NLDAS wind" &
      preflight_inventory$exists
  )
) {

  validation <- tryCatch({

    wind_brick <- raster::brick(
      preflight_inventory$resolved_file[i]
    )

    layer_names <- names(
      wind_brick
    )

    list(
      n_layers = raster::nlayers(
        wind_brick
      ),
      layer_names = paste(
        layer_names,
        collapse = "|"
      ),
      valid =
        raster::nlayers(
          wind_brick
        ) ==
        2 &&
        all(
          c(
            "Wind_E",
            "Wind_N"
          ) %in%
            layer_names
        )
    )

  }, error = function(
      e) {

    list(
      n_layers = NA_integer_,
      layer_names = paste0(
        "ERROR: ",
        conditionMessage(
          e
        )
      ),
      valid = FALSE
    )
  })

  preflight_inventory$wind_layers[i] <- validation$n_layers
  preflight_inventory$wind_layer_names[i] <- validation$layer_names
  preflight_inventory$wind_valid[i] <- validation$valid
}

preflight_inventory[
  input_type !=
    "Annual NLDAS wind",
  wind_valid := NA
]

fwrite(
  preflight_inventory,
  preflight_inventory_file
)

cat(
  "\n========================================\n",
  "STRICT PREFLIGHT INVENTORY\n",
  "========================================\n",
  sep = ""
)

print(
  preflight_inventory[, .(
    input_type,
    year,
    exists,
    resolution_status,
    candidate_count,
    wind_layers,
    wind_layer_names,
    resolved_file
  )]
)

bad_nonwind <- preflight_inventory[
  input_type !=
    "Annual NLDAS wind" &
    !exists
]

bad_wind <- preflight_inventory[
  input_type ==
    "Annual NLDAS wind" &
    (
      !exists |
        is.na(
          wind_valid
        ) |
        !wind_valid
    )
]

ambiguous_wind <- preflight_inventory[
  input_type ==
    "Annual NLDAS wind" &
    candidate_count >
    1L
]

if (
  nrow(
    bad_nonwind
  ) >
    0 ||
    nrow(
      bad_wind
    ) >
    0 ||
    nrow(
      ambiguous_wind
    ) >
    0
) {
  stop(
    "STRICT PREFLIGHT FAILED. No model years were processed.\n\n",
    "Review:\n",
    preflight_inventory_file,
    "\n\nMissing/invalid non-wind inputs:\n",
    if (
      nrow(
        bad_nonwind
      ) >
        0
    ) {
      paste(
        bad_nonwind$expected_file,
        collapse = "\n"
      )
    } else {
      "(none)"
    },
    "\n\nMissing/invalid/ambiguous wind inputs:\n",
    if (
      nrow(
        bad_wind
      ) >
        0
    ) {
      paste(
        paste0(
          bad_wind$year,
          ": ",
          bad_wind$resolution_status,
          " | ",
          bad_wind$resolved_file
        ),
        collapse = "\n"
      )
    } else {
      "(none)"
    },
    "\n\nA complete 2011-2020 result cannot be produced until every year passes."
  )
}

wind_files <- setNames(
  wind_inventory$resolved_file,
  as.character(
    wind_inventory$year
  )
)

if (
  !identical(
    as.integer(
      names(
        wind_files
      )
    ),
    as.integer(
      years
    )
  )
) {
  stop(
    "Internal error: resolved wind-file years do not match 2011-2020."
  )
}

cat(
  "\nSTRICT PREFLIGHT PASSED: all 2011-2020 inputs are present and readable.\n",
  sep = ""
)

# Remove only final products from a prior strict run after preflight succeeds.
unlink(
  c(
    pdf_file,
    png_file,
    summary_qc_file,
    year_qc_file,
    parent_qc_file,
    plot_data_file,
    gasoline_output_file,
    diesel_output_file,
    completion_file
  ),
  force = TRUE
)

# Display scales
pm25_min <- 0
pm25_max <- 1.0
pm25_breaks <- seq(
  0,
  1.0,
  by = 0.2
)

difference_upper_quantile <- 0.995

# Wind-speed modifier parameters retained from the earlier model.
U0_ms <- 0.1
beta_wind <- 0.5
angular_offset_rad <- 0.05

# A vector magnitude below this threshold is treated as directionally calm.
calm_vector_threshold_ms <- 0.05

# Numerical tolerances
small_value <- 1e-12
mapping_tolerance <- 1e-8

# Figure dimensions
fig_width_in <- 12.2
fig_height_in <- 6.75
png_dpi <- 600

# ------------------------------------------------------------------------------
# 2. Read Figure 4-3 templates
# ------------------------------------------------------------------------------

read_roadiness_template <- function(
    file_path,
    source_code,
    source_name) {

  dt <- read_fst(
    file_path,
    as.data.table = TRUE
  )

  required_columns <- c(
    "x",
    "y",
    "parent_id",
    "pm12",
    "roadiness",
    "weight",
    "pm1"
  )

  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )

  if (
    length(
      missing_columns
    ) > 0
  ) {
    stop(
      "Figure 4-3 template is missing required columns:\n",
      file_path,
      "\nMissing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  out <- dt[, .(
    x = as.numeric(x),
    y = as.numeric(y),
    parent_id = as.integer(parent_id),
    decade_pm12_template = as.numeric(pm12),
    roadiness = as.numeric(roadiness),
    roadiness_weight = as.numeric(weight),
    decade_roadiness_pm = as.numeric(pm1)
  )]

  out <- out[
    is.finite(x) &
      is.finite(y) &
      is.finite(parent_id) &
      is.finite(decade_pm12_template) &
      is.finite(roadiness) &
      is.finite(roadiness_weight) &
      is.finite(decade_roadiness_pm)
  ]

  if (
    nrow(
      out
    ) == 0
  ) {
    stop(
      "No usable records were found in:\n",
      file_path
    )
  }

  if (
    anyDuplicated(
      out[, .(
        x,
        y
      )]
    ) > 0
  ) {
    stop(
      "Duplicate 1-km coordinates were found in:\n",
      file_path
    )
  }

  out[, source_code := source_code]
  out[, source_name := source_name]

  out
}

gasoline_template <- read_roadiness_template(
  file_path = gasoline_template_file,
  source_code = "ONR",
  source_name = "Gasoline vehicles"
)

diesel_template <- read_roadiness_template(
  file_path = diesel_template_file,
  source_code = "NRD",
  source_name = "Diesel vehicles"
)

gasoline_xy <- gasoline_template[
  order(
    x,
    y
  ),
  .(
    x,
    y,
    parent_id
  )
]

diesel_xy <- diesel_template[
  order(
    x,
    y
  ),
  .(
    x,
    y,
    parent_id
  )
]

if (
  !identical(
    gasoline_xy,
    diesel_xy
  )
) {
  stop(
    "Gasoline and diesel Figure 4-3 templates do not use identical ",
    "1-km coordinates and parent assignments."
  )
}

template_xy <- gasoline_template[, .(
  x,
  y,
  parent_id
)]

# Check that roadiness-only weights average to one within each parent cell.
template_weight_qc <- rbindlist(
  list(
    gasoline_template[, .(
      mean_weight = mean(
        roadiness_weight,
        na.rm = TRUE
      )
    ), by = parent_id][, source_code := "ONR"],

    diesel_template[, .(
      mean_weight = mean(
        roadiness_weight,
        na.rm = TRUE
      )
    ), by = parent_id][, source_code := "NRD"]
  ),
  use.names = TRUE
)

maximum_template_weight_error <- max(
  abs(
    template_weight_qc$mean_weight -
      1
  ),
  na.rm = TRUE
)

if (
  maximum_template_weight_error >
    1e-7
) {
  warning(
    "Figure 4-3 roadiness-only weights do not average exactly to one ",
    "within all parent cells. Maximum error = ",
    signif(
      maximum_template_weight_error,
      6
    )
  )
}

# ------------------------------------------------------------------------------
# 3. Virginia boundary and 1-km point coordinates in WGS84
# ------------------------------------------------------------------------------

states_ll <- USAboundaries::us_states(
  resolution = "low"
)

virginia_ll <- states_ll[
  states_ll$state_abbr == "VA",
]

if (
  nrow(
    virginia_ll
  ) != 1
) {
  stop(
    "Virginia was not uniquely identified."
  )
}

virginia_lcc <- st_transform(
  virginia_ll,
  crs = st_crs(
    p4s
  )
)

template_points_lcc <- st_as_sf(
  template_xy,
  coords = c(
    "x",
    "y"
  ),
  crs = st_crs(
    p4s
  ),
  remove = FALSE
)

template_points_ll <- st_transform(
  template_points_lcc,
  4326
)

template_lonlat <- st_coordinates(
  template_points_ll
)

# ------------------------------------------------------------------------------
# 4. Reconstruct the same native parent grid used by Figure 4-3
# ------------------------------------------------------------------------------

read_decade_source <- function(
    file_path) {

  dt <- as.data.table(
    readRDS(
      file_path
    )
  )

  required_columns <- c(
    "x",
    "y",
    "value"
  )

  missing_columns <- setdiff(
    required_columns,
    names(
      dt
    )
  )

  if (
    length(
      missing_columns
    ) > 0
  ) {
    stop(
      "Decade CMAQ file is missing columns:\n",
      file_path,
      "\nMissing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  out <- dt[, .(
    lon = as.numeric(x),
    lat = as.numeric(y),
    concentration = as.numeric(value)
  )]

  out[
    is.finite(
      lon
    ) &
      is.finite(
        lat
      ) &
      is.finite(
        concentration
      )
  ]
}

gasoline_decade <- read_decade_source(
  gasoline_decade_file
)

diesel_decade <- read_decade_source(
  diesel_decade_file
)

va_bbox_ll <- st_bbox(
  virginia_ll
)

cmaq_margin_deg <- 1.5

cmaq_grid <- unique(
  gasoline_decade[
    lon >=
      as.numeric(
        va_bbox_ll["xmin"]
      ) -
      cmaq_margin_deg &
      lon <=
      as.numeric(
        va_bbox_ll["xmax"]
      ) +
      cmaq_margin_deg &
      lat >=
      as.numeric(
        va_bbox_ll["ymin"]
      ) -
      cmaq_margin_deg &
      lat <=
      as.numeric(
        va_bbox_ll["ymax"]
      ) +
      cmaq_margin_deg,
    .(
      lon,
      lat
    )
  ]
)

setorder(
  cmaq_grid,
  lon,
  lat
)

cmaq_grid[, parent_id := .I]

if (
  max(
    template_xy$parent_id,
    na.rm = TRUE
  ) >
    nrow(
      cmaq_grid
    )
) {
  stop(
    "Figure 4-3 parent IDs exceed the reconstructed CMAQ grid size."
  )
}

# Verify parent mapping against the stored Figure 4-3 decade concentrations.
verify_template_mapping <- function(
    decade_dt,
    template_dt,
    source_name) {

  parent_values <- merge(
    cmaq_grid,
    decade_dt[, .(
      lon,
      lat,
      expected_pm12 = concentration
    )],
    by = c(
      "lon",
      "lat"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  setorder(
    parent_values,
    parent_id
  )

  if (
    anyNA(
      parent_values$expected_pm12
    )
  ) {
    stop(
      "Missing decade parent concentrations while checking ",
      source_name,
      "."
    )
  }

  mapped <- parent_values$expected_pm12[
    template_dt$parent_id
  ]

  maximum_difference <- max(
    abs(
      mapped -
        template_dt$decade_pm12_template
    ),
    na.rm = TRUE
  )

  if (
    maximum_difference >
      mapping_tolerance
  ) {
    stop(
      "Reconstructed parent mapping does not reproduce the Figure 4-3 ",
      source_name,
      " template. Maximum difference = ",
      maximum_difference
    )
  }

  invisible(
    parent_values
  )
}

verify_template_mapping(
  decade_dt = gasoline_decade,
  template_dt = gasoline_template,
  source_name = "gasoline"
)

verify_template_mapping(
  decade_dt = diesel_decade,
  template_dt = diesel_template,
  source_name = "diesel"
)

# ------------------------------------------------------------------------------
# 5. Precompute eight-neighbor indices on the 1-km grid
# ------------------------------------------------------------------------------

detect_grid_spacing <- function(
    x) {

  unique_values <- sort(
    unique(
      x
    )
  )

  differences <- diff(
    unique_values
  )

  differences <- differences[
    is.finite(
      differences
    ) &
      differences >
      0
  ]

  median(
    differences,
    na.rm = TRUE
  )
}

dx_grid <- detect_grid_spacing(
  template_xy$x
)

dy_grid <- detect_grid_spacing(
  template_xy$y
)

if (
  !is.finite(
    dx_grid
  ) ||
    !is.finite(
      dy_grid
    )
) {
  stop(
    "Could not determine 1-km grid spacing."
  )
}

cat(
  "\nDetected grid spacing: ",
  dx_grid,
  " x ",
  dy_grid,
  " map units\n",
  sep = ""
)

coordinate_key <- function(
    x,
    y) {

  paste0(
    sprintf(
      "%.3f",
      x
    ),
    "|",
    sprintf(
      "%.3f",
      y
    )
  )
}

template_keys <- coordinate_key(
  template_xy$x,
  template_xy$y
)

directions <- data.table(
  name = c(
    "E",
    "NE",
    "N",
    "NW",
    "W",
    "SW",
    "S",
    "SE"
  ),
  sx = c(
    1,
    1,
    0,
    -1,
    -1,
    -1,
    0,
    1
  ),
  sy = c(
    0,
    1,
    1,
    1,
    0,
    -1,
    -1,
    -1
  ),
  angle = c(
    0,
    pi / 4,
    pi / 2,
    3 * pi / 4,
    pi,
    5 * pi / 4,
    3 * pi / 2,
    7 * pi / 4
  )
)

directions[, distance_km :=
             sqrt(
               sx^2 +
                 sy^2
             )]

neighbor_indices <- vector(
  "list",
  nrow(
    directions
  )
)

names(
  neighbor_indices
) <- directions$name

for (
  k in seq_len(
    nrow(
      directions
    )
  )
) {

  neighbor_key <- coordinate_key(
    template_xy$x +
      directions$sx[k] *
      dx_grid,
    template_xy$y +
      directions$sy[k] *
      dy_grid
  )

  neighbor_indices[[k]] <- match(
    neighbor_key,
    template_keys
  )
}

# ------------------------------------------------------------------------------
# 6. Read annual CMAQ parent values
# ------------------------------------------------------------------------------

read_annual_parent_values <- function(
    source_code,
    year) {

  file_path <- file.path(
    cmaq_dir,
    sprintf(
      "ANNUAL_%s_mean_%d.rds",
      source_code,
      year
    )
  )

  dt <- as.data.table(
    readRDS(
      file_path
    )
  )

  required_columns <- c(
    "x",
    "y",
    "value"
  )

  missing_columns <- setdiff(
    required_columns,
    names(
      dt
    )
  )

  if (
    length(
      missing_columns
    ) > 0
  ) {
    stop(
      "Annual CMAQ file is missing columns:\n",
      file_path,
      "\nMissing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  annual_dt <- dt[, .(
    lon = as.numeric(x),
    lat = as.numeric(y),
    pm12 = as.numeric(value)
  )]

  annual_dt <- annual_dt[
    is.finite(
      lon
    ) &
      is.finite(
        lat
      ) &
      is.finite(
        pm12
      )
  ]

  parent_values <- merge(
    cmaq_grid,
    annual_dt,
    by = c(
      "lon",
      "lat"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  setorder(
    parent_values,
    parent_id
  )

  if (
    anyNA(
      parent_values$pm12
    )
  ) {
    stop(
      "Some reconstructed parent cells are missing annual ",
      source_code,
      " concentrations for ",
      year,
      "."
    )
  }

  parent_values
}

# ------------------------------------------------------------------------------
# 7. Extract annual wind components
# ------------------------------------------------------------------------------

extract_annual_wind <- function(
    year) {

  wind_file <- unname(
    wind_files[
      as.character(
        year
      )
    ]
  )

  if (
    length(
      wind_file
    ) !=
      1 ||
      is.na(
        wind_file
      ) ||
      !file.exists(
        wind_file
      )
  ) {
    stop(
      "Resolved wind file is unavailable for ",
      year,
      ": ",
      wind_file
    )
  }

  wind_brick <- raster::brick(
    wind_file
  )

  if (
    raster::nlayers(
      wind_brick
    ) != 2
  ) {
    stop(
      "Expected two wind layers in:\n",
      wind_file
    )
  }

  layer_names <- names(
    wind_brick
  )

  if (
    !all(
      c(
        "Wind_E",
        "Wind_N"
      ) %in%
        layer_names
    )
  ) {
    stop(
      "Wind file does not contain Wind_E and Wind_N:\n",
      wind_file,
      "\nLayers: ",
      paste(
        layer_names,
        collapse = ", "
      )
    )
  }

  # Use nearest-cell extraction, matching the earlier implementation.
  u <- raster::extract(
    wind_brick[["Wind_E"]],
    template_lonlat,
    method = "simple"
  )

  v <- raster::extract(
    wind_brick[["Wind_N"]],
    template_lonlat,
    method = "simple"
  )

  if (
    any(
      !is.finite(
        u
      ) |
        !is.finite(
          v
        )
    )
  ) {
    stop(
      "Non-finite wind components remained after extraction for ",
      year,
      "."
    )
  }

  wind_vector_magnitude <- sqrt(
    u^2 +
      v^2
  )

  theta_to <- atan2(
    v,
    u
  )

  theta_upwind <- (
    theta_to +
      pi
  ) %%
    (
      2 *
        pi
    )

  data.table(
    x = template_xy$x,
    y = template_xy$y,
    u = u,
    v = v,
    vector_magnitude = wind_vector_magnitude,
    theta_to = theta_to,
    theta_upwind = theta_upwind,
    direction_valid =
      wind_vector_magnitude >=
      calm_vector_threshold_ms
  )
}

# ------------------------------------------------------------------------------
# 8. Wind direction + wind speed modification
# ------------------------------------------------------------------------------

calculate_wind_weight <- function(
    template_dt,
    wind_dt) {

  road_local <- as.numeric(
    template_dt$roadiness
  )

  n_cells <- nrow(
    template_dt
  )

  # Original formulation: all eight directions receive a positive angular
  # weight, strongest for the closest direction to the upwind bearing.
  numerator <- numeric(
    n_cells
  )

  denominator <- numeric(
    n_cells
  )

  for (
    k in seq_len(
      nrow(
        directions
      )
    )
  ) {

    neighbor_index <- neighbor_indices[[k]]

    # At the Virginia grid edge, the earlier code substituted the local
    # roadiness value when a neighboring cell was unavailable.
    neighbor_roadiness <- copy(
      road_local
    )

    valid_neighbor <- !is.na(
      neighbor_index
    )

    neighbor_roadiness[
      valid_neighbor
    ] <- road_local[
      neighbor_index[
        valid_neighbor
      ]
    ]

    angular_difference <- abs(
      atan2(
        sin(
          wind_dt$theta_upwind -
            directions$angle[k]
        ),
        cos(
          wind_dt$theta_upwind -
            directions$angle[k]
        )
      )
    )

    angular_weight <- 1 /
      (
        angular_difference +
          angular_offset_rad
      )

    usable <- is.finite(
      neighbor_roadiness
    ) &
      is.finite(
        angular_weight
      )

    numerator[
      usable
    ] <- numerator[
      usable
    ] +
      angular_weight[
        usable
      ] *
      neighbor_roadiness[
        usable
      ]

    denominator[
      usable
    ] <- denominator[
      usable
    ] +
      angular_weight[
        usable
      ]
  }

  upwind_roadiness <- numerator /
    (
      denominator +
        small_value
    )

  upwind_roadiness[
    !is.finite(
      upwind_roadiness
    )
  ] <- road_local[
    !is.finite(
      upwind_roadiness
    )
  ]

  domain_mean_wind_speed <- mean(
    wind_dt$vector_magnitude,
    na.rm = TRUE
  )

  if (
    !is.finite(
      domain_mean_wind_speed
    ) ||
      domain_mean_wind_speed <
      0
  ) {
    stop(
      "Invalid domain mean wind speed."
    )
  }

  speed_factor_beta <- (
    (
      domain_mean_wind_speed +
        small_value
    ) /
      (
        wind_dt$vector_magnitude +
          U0_ms
      )
  ) ^
    beta_wind

  speed_factor_beta[
    !is.finite(
      speed_factor_beta
    ) |
      speed_factor_beta <
      0
  ] <- 1

  # Original model structure:
  # direction-adjusted upwind roadiness multiplied by the speed factor.
  raw_wind_weight <- upwind_roadiness *
    speed_factor_beta

  raw_wind_weight[
    !is.finite(
      raw_wind_weight
    ) |
      raw_wind_weight <
      0
  ] <- 0

  result <- data.table(
    x = template_dt$x,
    y = template_dt$y,
    parent_id = template_dt$parent_id,
    roadiness = road_local,
    upwind_roadiness = upwind_roadiness,
    domain_mean_wind_speed = domain_mean_wind_speed,
    speed_factor_beta = speed_factor_beta,
    raw_wind_weight = raw_wind_weight,
    u = wind_dt$u,
    v = wind_dt$v,
    vector_magnitude = wind_dt$vector_magnitude,
    theta_to = wind_dt$theta_to,
    theta_upwind = wind_dt$theta_upwind,
    direction_valid = wind_dt$direction_valid
  )

  result[, parent_wind_mean :=
           mean(
             raw_wind_weight,
             na.rm = TRUE
           ),
         by = parent_id]

  result[, wind_weight :=
           fifelse(
             is.finite(
               parent_wind_mean
             ) &
               parent_wind_mean >
               0,
             (
               raw_wind_weight +
                 small_value
             ) /
               (
                 parent_wind_mean +
                   small_value
               ),
             1
           )]

  if (
    any(
      !is.finite(
        result$wind_weight
      )
    )
  ) {
    stop(
      "Non-finite direction-and-speed wind weights were produced."
    )
  }

  result
}

# ------------------------------------------------------------------------------
# 9. Annual processing
# ------------------------------------------------------------------------------

annual_results <- list()
year_qc_list <- list()
parent_qc_list <- list()

for (
  year in years
) {

  cat(
    "\n========================================\n",
    "PROCESSING YEAR ",
    year,
    "\n",
    "========================================\n",
    sep = ""
  )

  wind_dt <- extract_annual_wind(
    year
  )

  # Domain-level vector summary.
  mean_u <- mean(
    wind_dt$u,
    na.rm = TRUE
  )

  mean_v <- mean(
    wind_dt$v,
    na.rm = TRUE
  )

  mean_vector_magnitude <- sqrt(
    mean_u^2 +
      mean_v^2
  )

  mean_to_bearing_deg <- (
    atan2(
      mean_u,
      mean_v
    ) *
      180 /
      pi +
      360
  ) %%
    360

  mean_from_bearing_deg <- (
    mean_to_bearing_deg +
      180
  ) %%
    360

  for (
    source_code in c(
      "ONR",
      "NRD"
    )
  ) {

    source_name <- if (
      source_code ==
        "ONR"
    ) {
      "Gasoline vehicles"
    } else {
      "Diesel vehicles"
    }

    template_dt <- if (
      source_code ==
        "ONR"
    ) {
      gasoline_template
    } else {
      diesel_template
    }

    parent_values <- read_annual_parent_values(
      source_code = source_code,
      year = year
    )

    wind_weight_dt <- calculate_wind_weight(
      template_dt = template_dt,
      wind_dt = wind_dt
    )

    annual_dt <- copy(
      wind_weight_dt
    )

    annual_dt[, pm12 :=
                parent_values$pm12[
                  parent_id
                ]]

    annual_dt[, roadiness_weight :=
                template_dt$roadiness_weight]

    annual_dt[, pm_roadiness :=
                pm12 *
                roadiness_weight]

    annual_dt[, pm_wind :=
                pm12 *
                wind_weight]

    annual_dt[, delta_wind :=
                pm_wind -
                pm_roadiness]

    annual_dt[, year := year]
    annual_dt[, source_code := source_code]
    annual_dt[, source_name := source_name]

    parent_check <- annual_dt[, .(
      pm12 = first(
        pm12
      ),
      mean_roadiness = mean(
        pm_roadiness,
        na.rm = TRUE
      ),
      mean_wind = mean(
        pm_wind,
        na.rm = TRUE
      ),
      roadiness_error =
        mean(
          pm_roadiness,
          na.rm = TRUE
        ) -
        first(
          pm12
        ),
      wind_error =
        mean(
          pm_wind,
          na.rm = TRUE
        ) -
        first(
          pm12
        ),
      n_1km_cells = .N
    ), by = parent_id]

    parent_check[, year := year]
    parent_check[, source_code := source_code]
    parent_check[, source_name := source_name]

    parent_qc_list[[
      paste0(
        source_code,
        "_",
        year
      )
    ]] <- parent_check

    year_qc <- data.table(
      year = year,
      source_code = source_code,
      source_name = source_name,
      n_1km_cells = nrow(
        annual_dt
      ),
      n_parent_cells = uniqueN(
        annual_dt$parent_id
      ),
      mean_pm_roadiness = mean(
        annual_dt$pm_roadiness,
        na.rm = TRUE
      ),
      mean_pm_wind = mean(
        annual_dt$pm_wind,
        na.rm = TRUE
      ),
      mean_delta = mean(
        annual_dt$delta_wind,
        na.rm = TRUE
      ),
      mean_absolute_delta = mean(
        abs(
          annual_dt$delta_wind
        ),
        na.rm = TRUE
      ),
      p95_absolute_delta = as.numeric(
        quantile(
          abs(
            annual_dt$delta_wind
          ),
          0.95,
          na.rm = TRUE
        )
      ),
      delta_p05 = as.numeric(
        quantile(
          annual_dt$delta_wind,
          0.05,
          na.rm = TRUE
        )
      ),
      delta_p95 = as.numeric(
        quantile(
          annual_dt$delta_wind,
          0.95,
          na.rm = TRUE
        )
      ),
      maximum_absolute_roadiness_conservation_error =
        max(
          abs(
            parent_check$roadiness_error
          ),
          na.rm = TRUE
        ),
      maximum_absolute_wind_conservation_error =
        max(
          abs(
            parent_check$wind_error
          ),
          na.rm = TRUE
        ),
      mean_wind_eastward_component = mean_u,
      mean_wind_northward_component = mean_v,
      domain_mean_vector_magnitude = mean_vector_magnitude,
      mean_local_wind_speed = mean(
        annual_dt$vector_magnitude,
        na.rm = TRUE
      ),
      mean_speed_factor_beta = mean(
        annual_dt$speed_factor_beta,
        na.rm = TRUE
      ),
      speed_factor_beta_p05 = as.numeric(
        quantile(
          annual_dt$speed_factor_beta,
          0.05,
          na.rm = TRUE
        )
      ),
      speed_factor_beta_p95 = as.numeric(
        quantile(
          annual_dt$speed_factor_beta,
          0.95,
          na.rm = TRUE
        )
      ),
      domain_mean_wind_to_bearing_deg = mean_to_bearing_deg,
      domain_mean_wind_from_bearing_deg = mean_from_bearing_deg,
      percent_directionally_calm =
        100 *
        mean(
          !wind_dt$direction_valid,
          na.rm = TRUE
        )
    )

    year_qc_list[[
      paste0(
        source_code,
        "_",
        year
      )
    ]] <- year_qc

    annual_results[[
      paste0(
        source_code,
        "_",
        year
      )
    ]] <- annual_dt[, .(
      x,
      y,
      parent_id,
      year,
      source_code,
      source_name,
      pm12,
      pm_roadiness,
      pm_wind,
      delta_wind,
      roadiness,
      upwind_roadiness,
      speed_factor_beta,
      raw_wind_weight,
      roadiness_weight,
      wind_weight,
      u,
      v,
      vector_magnitude,
      theta_to,
      theta_upwind
    )]
  }

  rm(
    wind_dt
  )

  gc()
}

annual_stack <- rbindlist(
  annual_results,
  use.names = TRUE
)

year_qc <- rbindlist(
  year_qc_list,
  use.names = TRUE,
  fill = TRUE
)

parent_qc <- rbindlist(
  parent_qc_list,
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------------------------
# 9A. Strict post-processing completeness check
# ------------------------------------------------------------------------------

completed_years <- sort(
  unique(
    annual_stack$year
  )
)

if (
  !identical(
    as.integer(
      completed_years
    ),
    as.integer(
      years
    )
  )
) {
  stop(
    "STRICT COMPLETENESS CHECK FAILED before any final output was written.\n",
    "Expected years: ",
    paste(
      years,
      collapse = ", "
    ),
    "\nCompleted years: ",
    paste(
      completed_years,
      collapse = ", "
    )
  )
}

source_year_counts <- annual_stack[, .(
  n_completed_years = uniqueN(
    year
  )
), by = source_code]

if (
  nrow(
    source_year_counts
  ) !=
    2 ||
    any(
      source_year_counts$n_completed_years !=
        length(
          years
        )
    )
) {
  stop(
    "STRICT COMPLETENESS CHECK FAILED: each source must contain exactly ",
    length(
      years
    ),
    " years.\n",
    paste(
      capture.output(
        print(
          source_year_counts
        )
      ),
      collapse = "\n"
    )
  )
}

expected_source_years <- CJ(
  source_code = c(
    "NRD",
    "ONR"
  ),
  year = years
)

actual_source_years <- unique(
  annual_stack[, .(
    source_code,
    year
  )]
)

setkey(
  expected_source_years,
  source_code,
  year
)

setkey(
  actual_source_years,
  source_code,
  year
)

missing_source_years <- fsetdiff(
  expected_source_years,
  actual_source_years
)

if (
  nrow(
    missing_source_years
  ) >
    0
) {
  stop(
    "STRICT COMPLETENESS CHECK FAILED. Missing source-year combinations:\n",
    paste(
      capture.output(
        print(
          missing_source_years
        )
      ),
      collapse = "\n"
    )
  )
}

cat(
  "\nSTRICT COMPLETENESS CHECK PASSED: both sources contain all 10 years.\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 10. 2011-2020 mean outputs and summary QC
# ------------------------------------------------------------------------------

decade_data <- annual_stack[, .(
  pm12 = mean(
    pm12,
    na.rm = TRUE
  ),
  pm_roadiness = mean(
    pm_roadiness,
    na.rm = TRUE
  ),
  pm_wind = mean(
    pm_wind,
    na.rm = TRUE
  ),
  delta_wind = mean(
    delta_wind,
    na.rm = TRUE
  ),
  roadiness = first(
    roadiness
  ),
  mean_upwind_roadiness = mean(
    upwind_roadiness,
    na.rm = TRUE
  ),
  mean_speed_factor_beta = mean(
    speed_factor_beta,
    na.rm = TRUE
  ),
  mean_raw_wind_weight = mean(
    raw_wind_weight,
    na.rm = TRUE
  ),
  roadiness_weight = first(
    roadiness_weight
  ),
  mean_wind_weight = mean(
    wind_weight,
    na.rm = TRUE
  ),
  mean_u = mean(
    u,
    na.rm = TRUE
  ),
  mean_v = mean(
    v,
    na.rm = TRUE
  ),
  mean_vector_magnitude = mean(
    vector_magnitude,
    na.rm = TRUE
  ),
  n_years = uniqueN(
    year
  )
), by = .(
  x,
  y,
  parent_id,
  source_code,
  source_name
)]

summary_qc <- decade_data[, .(
  n_1km_cells = .N,
  n_parent_cells = uniqueN(
    parent_id
  ),
  mean_roadiness = mean(
    pm_roadiness,
    na.rm = TRUE
  ),
  mean_wind = mean(
    pm_wind,
    na.rm = TRUE
  ),
  mean_delta = mean(
    delta_wind,
    na.rm = TRUE
  ),
  median_absolute_delta = median(
    abs(
      delta_wind
    ),
    na.rm = TRUE
  ),
  mean_absolute_delta = mean(
    abs(
      delta_wind
    ),
    na.rm = TRUE
  ),
  p95_absolute_delta = as.numeric(
    quantile(
      abs(
        delta_wind
      ),
      0.95,
      na.rm = TRUE
    )
  ),
  delta_p05 = as.numeric(
    quantile(
      delta_wind,
      0.05,
      na.rm = TRUE
    )
  ),
  delta_p95 = as.numeric(
    quantile(
      delta_wind,
      0.95,
      na.rm = TRUE
    )
  ),
  delta_minimum = min(
    delta_wind,
    na.rm = TRUE
  ),
  delta_maximum = max(
    delta_wind,
    na.rm = TRUE
  ),
  percent_absolute_change_above_0_01 =
    100 *
    mean(
      abs(
        delta_wind
      ) >=
        0.01,
      na.rm = TRUE
    ),
  percent_absolute_change_above_0_05 =
    100 *
    mean(
      abs(
        delta_wind
      ) >=
        0.05,
      na.rm = TRUE
    ),
  mean_speed_factor_beta = mean(
    mean_speed_factor_beta,
    na.rm = TRUE
  ),
  speed_factor_beta_p05 = as.numeric(
    quantile(
      mean_speed_factor_beta,
      0.05,
      na.rm = TRUE
    )
  ),
  speed_factor_beta_p95 = as.numeric(
    quantile(
      mean_speed_factor_beta,
      0.95,
      na.rm = TRUE
    )
  ),
  maximum_absolute_wind_conservation_error =
    max(
      abs(
        parent_qc[
          source_code ==
            .BY$source_code,
          wind_error
        ]
      ),
      na.rm = TRUE
    ),
  minimum_years = min(
    n_years,
    na.rm = TRUE
  ),
  maximum_years = max(
    n_years,
    na.rm = TRUE
  )
), by = .(
  source_code,
  source_name
)]

fwrite(
  summary_qc,
  summary_qc_file
)

fwrite(
  year_qc,
  year_qc_file
)

fwrite(
  parent_qc,
  parent_qc_file
)

saveRDS(
  list(
    decade_data = decade_data,
    annual_data = annual_stack,
    summary_qc = summary_qc,
    year_qc = year_qc,
    parent_qc = parent_qc,
    method = list(
      wind_components =
        "Annual signed Wind_E and Wind_N",
      direction_weight =
        "1 / (absolute angular difference + 0.05), all eight neighbors",
      raw_wind_modifier =
        "upwind roadiness multiplied by relative wind-speed factor",
      speed_factor =
        "((domain mean vector magnitude + eps) / (local vector magnitude + 0.1))^0.5",
      parent_normalization =
        "raw modifier divided by parent-cell arithmetic mean"
    )
  ),
  plot_data_file
)

write_fst(
  decade_data[
    source_code ==
      "ONR"
  ],
  gasoline_output_file
)

write_fst(
  decade_data[
    source_code ==
      "NRD"
  ],
  diesel_output_file
)

cat(
  "\n========================================\n",
  "FIGURE 4-4 QC SUMMARY\n",
  "========================================\n",
  sep = ""
)

print(
  summary_qc
)

cat(
  "\nAnnual wind summaries:\n",
  sep = ""
)

print(
  unique(
    year_qc[, .(
      year,
      mean_wind_eastward_component,
      mean_wind_northward_component,
      domain_mean_vector_magnitude,
      domain_mean_wind_to_bearing_deg,
      domain_mean_wind_from_bearing_deg,
      percent_directionally_calm
    )]
  )
)

# ------------------------------------------------------------------------------
# 11. Figure scales and Virginia extent
# ------------------------------------------------------------------------------

difference_limit <- as.numeric(
  quantile(
    abs(
      decade_data$delta_wind
    ),
    difference_upper_quantile,
    na.rm = TRUE
  )
)

if (
  !is.finite(
    difference_limit
  ) ||
    difference_limit <=
    0
) {
  difference_limit <- max(
    abs(
      decade_data$delta_wind
    ),
    na.rm = TRUE
  )
}

difference_breaks <- pretty(
  c(
    -difference_limit,
    difference_limit
  ),
  n = 5
)

difference_breaks <- difference_breaks[
  difference_breaks >=
    -difference_limit &
    difference_breaks <=
    difference_limit
]

if (
  length(
    difference_breaks
  ) <
    3
) {
  difference_breaks <- seq(
    -difference_limit,
    difference_limit,
    length.out = 5
  )
}

va_bbox_lcc <- st_bbox(
  virginia_lcc
)

map_padding_m <- 10000

va_xlim <- c(
  as.numeric(
    va_bbox_lcc["xmin"]
  ) -
    map_padding_m,
  as.numeric(
    va_bbox_lcc["xmax"]
  ) +
    map_padding_m
)

va_ylim <- c(
  as.numeric(
    va_bbox_lcc["ymin"]
  ) -
    map_padding_m,
  as.numeric(
    va_bbox_lcc["ymax"]
  ) +
    map_padding_m
)

cat(
  "\nAbsolute PM2.5 display range: 0-1.0 ug/m3\n",
  "Difference display range: +/-",
  signif(
    difference_limit,
    5
  ),
  " ug/m3\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 12. Publication-style map panels
# ------------------------------------------------------------------------------

make_absolute_panel <- function(
    data,
    value_column,
    show_legend = FALSE) {

  panel_dt <- copy(
    data
  )

  panel_dt[, value_to_plot :=
             as.numeric(
               get(
                 value_column
               )
             )]

  ggplot() +

    geom_raster(
      data = panel_dt,
      aes(
        x = x,
        y = y,
        fill = value_to_plot
      ),
      interpolate = FALSE
    ) +

    geom_sf(
      data = virginia_lcc,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey25",
      linewidth = 0.34
    ) +

    coord_sf(
      crs = st_crs(
        p4s
      ),
      xlim = va_xlim,
      ylim = va_ylim,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +

    scale_fill_viridis_c(
      option = "E",
      direction = 1,
      limits = c(
        pm25_min,
        pm25_max
      ),
      breaks = pm25_breaks,
      labels = scales::label_number(
        accuracy = 0.1
      ),
      oob = scales::squish,
      na.value = "white",
      name = expression(
        PM[2.5] ~
          (mu * g ~ m^{-3})
      ),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = grid::unit(
          80,
          "mm"
        ),
        barheight = grid::unit(
          4.0,
          "mm"
        ),
        frame.colour = "grey35",
        frame.linewidth = 0.30,
        ticks.colour = "grey35"
      )
    ) +

    labs(
      x = NULL,
      y = NULL
    ) +

    theme_void(
      base_size = 10,
      base_family = "sans"
    ) +

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
        color = "grey38",
        linewidth = 0.30
      ),

      legend.position = if (
        show_legend
      ) {
        "bottom"
      } else {
        "none"
      },

      legend.justification = "center",

      legend.title = element_text(
        size = 9.1
      ),

      legend.text = element_text(
        size = 8.3
      ),

      legend.margin = margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),

      legend.box.margin = margin(
        0,
        0,
        0,
        0
      ),

      plot.margin = margin(
        1,
        4,
        3,
        4,
        unit = "pt"
      )
    )
}

make_difference_panel <- function(
    data,
    show_legend = FALSE) {

  ggplot() +

    geom_raster(
      data = data,
      aes(
        x = x,
        y = y,
        fill = delta_wind
      ),
      interpolate = FALSE
    ) +

    geom_sf(
      data = virginia_lcc,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey25",
      linewidth = 0.34
    ) +

    coord_sf(
      crs = st_crs(
        p4s
      ),
      xlim = va_xlim,
      ylim = va_ylim,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +

    scale_fill_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -difference_limit,
        difference_limit
      ),
      breaks = difference_breaks,
      labels = scales::label_number(
        accuracy = if (
          difference_limit <
            0.1
        ) {
          0.01
        } else {
          0.1
        }
      ),
      oob = scales::squish,
      na.value = "white",
      name = expression(
        Delta * PM[2.5] ~
          (mu * g ~ m^{-3})
      ),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = grid::unit(
          46,
          "mm"
        ),
        barheight = grid::unit(
          4.0,
          "mm"
        ),
        frame.colour = "grey35",
        frame.linewidth = 0.30,
        ticks.colour = "grey35"
      )
    ) +

    labs(
      x = NULL,
      y = NULL
    ) +

    theme_void(
      base_size = 10,
      base_family = "sans"
    ) +

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
        color = "grey38",
        linewidth = 0.30
      ),

      legend.position = if (
        show_legend
      ) {
        "bottom"
      } else {
        "none"
      },

      legend.justification = "center",

      legend.title = element_text(
        size = 9.1
      ),

      legend.text = element_text(
        size = 8.3
      ),

      legend.margin = margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0
      ),

      legend.box.margin = margin(
        0,
        0,
        0,
        0
      ),

      plot.margin = margin(
        1,
        4,
        3,
        4,
        unit = "pt"
      )
    )
}

gasoline_plot_data <- decade_data[
  source_code ==
    "ONR"
]

diesel_plot_data <- decade_data[
  source_code ==
    "NRD"
]

p_a <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm_roadiness",
  show_legend = FALSE
)

p_b <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm_wind",
  show_legend = FALSE
)

p_c <- make_difference_panel(
  data = gasoline_plot_data,
  show_legend = FALSE
)

p_d <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm_roadiness",
  show_legend = FALSE
)

p_e <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm_wind",
  show_legend = FALSE
)

p_f <- make_difference_panel(
  data = diesel_plot_data,
  show_legend = FALSE
)

p_absolute_legend <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm_wind",
  show_legend = TRUE
) +
  theme(
    panel.border = element_blank()
  )

p_difference_legend <- make_difference_panel(
  data = gasoline_plot_data,
  show_legend = TRUE
) +
  theme(
    panel.border = element_blank()
  )

get_legend <- function(
    plot_object) {

  plot_grob <- ggplotGrob(
    plot_object
  )

  guide_ids <- grep(
    "^guide-box",
    plot_grob$layout$name
  )

  if (
    length(
      guide_ids
    ) ==
      0
  ) {
    stop(
      "No legend guide-box was found."
    )
  }

  for (
    i in guide_ids
  ) {

    candidate <- plot_grob$grobs[[i]]

    if (
      !inherits(
        candidate,
        "zeroGrob"
      )
    ) {
      return(
        candidate
      )
    }
  }

  stop(
    "All legend guide-box objects were empty."
  )
}

absolute_legend_grob <- get_legend(
  p_absolute_legend
)

difference_legend_grob <- get_legend(
  p_difference_legend
)

g_a <- ggplotGrob(
  p_a
)

g_b <- ggplotGrob(
  p_b
)

g_c <- ggplotGrob(
  p_c
)

g_d <- ggplotGrob(
  p_d
)

g_e <- ggplotGrob(
  p_e
)

g_f <- ggplotGrob(
  p_f
)

# ------------------------------------------------------------------------------
# 13. Composite figure drawing
# ------------------------------------------------------------------------------

draw_figure <- function() {

  grid::grid.newpage()

  outer_layout <- grid::grid.layout(
    nrow = 3,
    ncol = 3,
    widths = grid::unit(
      c(
        1,
        1,
        1
      ),
      "null"
    ),
    heights = grid::unit.c(
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        0.76,
        "in"
      )
    )
  )

  grid::pushViewport(
    grid::viewport(
      layout = outer_layout
    )
  )

  draw_panel_with_tag <- function(
      grob_object,
      tag,
      row,
      column) {

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = row,
        layout.pos.col = column
      )
    )

    inner_layout <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(
          0.17,
          "in"
        ),
        grid::unit(
          1,
          "null"
        )
      )
    )

    grid::pushViewport(
      grid::viewport(
        layout = inner_layout
      )
    )

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    grid::grid.text(
      label = tag,
      x = grid::unit(
        1.5,
        "mm"
      ),
      y = grid::unit(
        0.5,
        "npc"
      ),
      just = c(
        "left",
        "center"
      ),
      gp = grid::gpar(
        fontfamily = "sans",
        fontface = "bold",
        fontsize = 10.5
      )
    )

    grid::popViewport()

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
      )
    )

    grid::grid.draw(
      grob_object
    )

    grid::popViewport(
      3
    )
  }

  draw_panel_with_tag(
    g_a,
    "(a)",
    1,
    1
  )

  draw_panel_with_tag(
    g_b,
    "(b)",
    1,
    2
  )

  draw_panel_with_tag(
    g_c,
    "(c)",
    1,
    3
  )

  draw_panel_with_tag(
    g_d,
    "(d)",
    2,
    1
  )

  draw_panel_with_tag(
    g_e,
    "(e)",
    2,
    2
  )

  draw_panel_with_tag(
    g_f,
    "(f)",
    2,
    3
  )

  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 3,
      layout.pos.col = 1:2
    )
  )

  grid::grid.draw(
    absolute_legend_grob
  )

  grid::popViewport()

  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 3,
      layout.pos.col = 3
    )
  )

  grid::grid.draw(
    difference_legend_grob
  )

  grid::popViewport()
  grid::popViewport()
}

# ------------------------------------------------------------------------------
# 14. Save vector PDF
# ------------------------------------------------------------------------------

unlink(
  c(
    pdf_file,
    png_file
  ),
  force = TRUE
)

grDevices::pdf(
  file = pdf_file,
  width = fig_width_in,
  height = fig_height_in,
  onefile = FALSE,
  useDingbats = FALSE,
  paper = "special",
  bg = "white",
  family = "Helvetica"
)

draw_figure()

invisible(
  grDevices::dev.off()
)

Sys.sleep(
  1
)

if (
  !file.exists(
    pdf_file
  ) ||
    !is.finite(
      file.info(
        pdf_file
      )$size
    ) ||
    file.info(
      pdf_file
    )$size <=
    1000
) {
  stop(
    "PDF creation failed:\n",
    pdf_file
  )
}

cat(
  "\nVector PDF created:\n",
  pdf_file,
  "\nSize: ",
  round(
    file.info(
      pdf_file
    )$size /
      1024^2,
    3
  ),
  " MB\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 15. Hopper-safe PNG creation with fallbacks
# ------------------------------------------------------------------------------

png_is_valid <- function(
    path) {

  file.exists(
    path
  ) &&
    is.finite(
      file.info(
        path
      )$size
    ) &&
    file.info(
      path
    )$size >
    1000
}

close_extra_devices <- function(
    start_device) {

  while (
    grDevices::dev.cur() >
      1 &&
      grDevices::dev.cur() !=
      start_device
  ) {
    try(
      grDevices::dev.off(),
      silent = TRUE
    )
  }
}

try_png_device <- function(
    method_name,
    open_device) {

  if (
    file.exists(
      png_file
    )
  ) {
    unlink(
      png_file,
      force = TRUE
    )
  }

  start_device <- grDevices::dev.cur()

  result <- tryCatch({

    open_device()

    draw_figure()

    invisible(
      grDevices::dev.off()
    )

    Sys.sleep(
      1
    )

    png_is_valid(
      png_file
    )

  }, error = function(e) {

    close_extra_devices(
      start_device
    )

    message(
      method_name,
      " failed: ",
      conditionMessage(
        e
      )
    )

    FALSE
  })

  if (
    isTRUE(
      result
    )
  ) {
    message(
      "PNG method succeeded: ",
      method_name
    )
  }

  isTRUE(
    result
  )
}

png_method <- NA_character_

# Method 1: base-R cairo-png
if (
  isTRUE(
    capabilities(
      "cairo"
    )
  ) &&
    is.na(
      png_method
    )
) {

  ok <- try_png_device(
    "base grDevices::png(type='cairo-png')",
    function() {

      grDevices::png(
        filename = png_file,
        width = as.integer(
          round(
            fig_width_in *
              png_dpi
          )
        ),
        height = as.integer(
          round(
            fig_height_in *
              png_dpi
          )
        ),
        units = "px",
        res = png_dpi,
        type = "cairo-png",
        bg = "white",
        pointsize = 12
      )
    }
  )

  if (
    ok
  ) {
    png_method <- "base cairo-png"
  }
}

# Method 2: ragg
if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  ) &&
    is.na(
      png_method
    )
) {

  ok <- try_png_device(
    "ragg::agg_png",
    function() {

      ragg::agg_png(
        filename = png_file,
        width = fig_width_in,
        height = fig_height_in,
        units = "in",
        res = png_dpi,
        background = "white"
      )
    }
  )

  if (
    ok
  ) {
    png_method <- "ragg"
  }
}

# Method 3: Cairo package
if (
  requireNamespace(
    "Cairo",
    quietly = TRUE
  ) &&
    is.na(
      png_method
    )
) {

  ok <- try_png_device(
    "Cairo::CairoPNG",
    function() {

      Cairo::CairoPNG(
        filename = png_file,
        width = as.integer(
          round(
            fig_width_in *
              png_dpi
          )
        ),
        height = as.integer(
          round(
            fig_height_in *
              png_dpi
          )
        ),
        unit = "px",
        pointsize = 12,
        bg = "white",
        res = png_dpi
      )
    }
  )

  if (
    ok
  ) {
    png_method <- "Cairo package"
  }
}

# Method 4: base bitmap / Ghostscript
if (
  is.na(
    png_method
  )
) {

  ok <- try_png_device(
    "grDevices::bitmap(type='png16m')",
    function() {

      grDevices::bitmap(
        file = png_file,
        type = "png16m",
        width = fig_width_in,
        height = fig_height_in,
        res = png_dpi,
        pointsize = 12
      )
    }
  )

  if (
    ok
  ) {
    png_method <- "base bitmap / Ghostscript"
  }
}

run_external <- function(
    command,
    args,
    method_name) {

  if (
    file.exists(
      png_file
    )
  ) {
    unlink(
      png_file,
      force = TRUE
    )
  }

  output <- tryCatch(
    system2(
      command = command,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      structure(
        conditionMessage(
          e
        ),
        status = 999L
      )
    }
  )

  status <- attr(
    output,
    "status"
  )

  if (
    is.null(
      status
    )
  ) {
    status <- 0L
  }

  if (
    status !=
      0L
  ) {
    message(
      method_name,
      " failed with status ",
      status,
      ":\n",
      paste(
        output,
        collapse = "\n"
      )
    )
  }

  Sys.sleep(
    1
  )

  if (
    status ==
      0L &&
      png_is_valid(
        png_file
      )
  ) {
    message(
      "PNG method succeeded: ",
      method_name
    )

    return(
      TRUE
    )
  }

  FALSE
}

find_existing_command <- function(
    command_name,
    fixed_paths = character()) {

  candidates <- unique(
    c(
      unname(
        Sys.which(
          command_name
        )
      ),
      fixed_paths
    )
  )

  candidates <- candidates[
    nzchar(
      candidates
    )
  ]

  candidates[
    file.exists(
      candidates
    )
  ][1]
}

# Method 5: Ghostscript
if (
  is.na(
    png_method
  )
) {

  gs_bin <- find_existing_command(
    "gs",
    c(
      "/usr/bin/gs",
      "/bin/gs",
      Sys.getenv(
        "R_GSCMD"
      )
    )
  )

  if (
    length(
      gs_bin
    ) ==
      1 &&
      !is.na(
        gs_bin
      )
  ) {

    ok <- run_external(
      command = gs_bin,
      args = c(
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=png16m",
        paste0(
          "-r",
          png_dpi
        ),
        "-dTextAlphaBits=4",
        "-dGraphicsAlphaBits=4",
        paste0(
          "-sOutputFile=",
          png_file
        ),
        pdf_file
      ),
      method_name = "Ghostscript"
    )

    if (
      ok
    ) {
      png_method <- "Ghostscript"
    }
  }
}

# Method 6: pdftocairo
if (
  is.na(
    png_method
  )
) {

  pdftocairo_bin <- find_existing_command(
    "pdftocairo",
    c(
      "/usr/bin/pdftocairo",
      "/bin/pdftocairo"
    )
  )

  if (
    length(
      pdftocairo_bin
    ) ==
      1 &&
      !is.na(
        pdftocairo_bin
      )
  ) {

    output_prefix <- tools::file_path_sans_ext(
      png_file
    )

    ok <- run_external(
      command = pdftocairo_bin,
      args = c(
        "-png",
        "-singlefile",
        "-r",
        as.character(
          png_dpi
        ),
        pdf_file,
        output_prefix
      ),
      method_name = "pdftocairo"
    )

    if (
      ok
    ) {
      png_method <- "pdftocairo"
    }
  }
}

# Method 7: pdftoppm
if (
  is.na(
    png_method
  )
) {

  pdftoppm_bin <- find_existing_command(
    "pdftoppm",
    c(
      "/usr/bin/pdftoppm",
      "/bin/pdftoppm"
    )
  )

  if (
    length(
      pdftoppm_bin
    ) ==
      1 &&
      !is.na(
        pdftoppm_bin
      )
  ) {

    output_prefix <- tools::file_path_sans_ext(
      png_file
    )

    ok <- run_external(
      command = pdftoppm_bin,
      args = c(
        "-png",
        "-r",
        as.character(
          png_dpi
        ),
        "-singlefile",
        pdf_file,
        output_prefix
      ),
      method_name = "pdftoppm"
    )

    if (
      ok
    ) {
      png_method <- "pdftoppm"
    }
  }
}

# Method 8: Python PyMuPDF (fitz), if available.
if (
  is.na(
    png_method
  )
) {

  python_bin <- unname(
    Sys.which(
      "python3"
    )
  )

  if (
    nzchar(
      python_bin
    ) &&
      file.exists(
        python_bin
      )
  ) {

    python_script <- file.path(
      out_dir,
      "_render_figure_4_4_with_pymupdf.py"
    )

    writeLines(
      c(
        "import sys",
        "import fitz",
        "pdf_file = sys.argv[1]",
        "png_file = sys.argv[2]",
        "dpi = float(sys.argv[3])",
        "doc = fitz.open(pdf_file)",
        "page = doc.load_page(0)",
        "matrix = fitz.Matrix(dpi / 72.0, dpi / 72.0)",
        "pix = page.get_pixmap(matrix=matrix, alpha=False)",
        "pix.save(png_file)",
        "doc.close()"
      ),
      con = python_script
    )

    ok <- run_external(
      command = python_bin,
      args = c(
        python_script,
        pdf_file,
        png_file,
        as.character(
          png_dpi
        )
      ),
      method_name = "Python PyMuPDF"
    )

    unlink(
      python_script,
      force = TRUE
    )

    if (
      ok
    ) {
      png_method <- "Python PyMuPDF"
    }
  }
}

# Method 9: ImageMagick
if (
  is.na(
    png_method
  )
) {

  magick_bin <- find_existing_command(
    "magick",
    c(
      "/usr/bin/magick",
      "/bin/magick"
    )
  )

  if (
    length(
      magick_bin
    ) ==
      1 &&
      !is.na(
        magick_bin
      )
  ) {

    ok <- run_external(
      command = magick_bin,
      args = c(
        "-density",
        as.character(
          png_dpi
        ),
        paste0(
          pdf_file,
          "[0]"
        ),
        "-background",
        "white",
        "-alpha",
        "remove",
        "-alpha",
        "off",
        png_file
      ),
      method_name = "ImageMagick magick"
    )

    if (
      ok
    ) {
      png_method <- "ImageMagick magick"
    }
  }
}

if (
  is.na(
    png_method
  )
) {

  convert_bin <- find_existing_command(
    "convert",
    c(
      "/usr/bin/convert",
      "/bin/convert"
    )
  )

  if (
    length(
      convert_bin
    ) ==
      1 &&
      !is.na(
        convert_bin
      )
  ) {

    ok <- run_external(
      command = convert_bin,
      args = c(
        "-density",
        as.character(
          png_dpi
        ),
        paste0(
          pdf_file,
          "[0]"
        ),
        "-background",
        "white",
        "-alpha",
        "remove",
        "-alpha",
        "off",
        png_file
      ),
      method_name = "ImageMagick convert"
    )

    if (
      ok
    ) {
      png_method <- "ImageMagick convert"
    }
  }
}

# Final fallback through Hopper module/login shell.
if (
  is.na(
    png_method
  ) &&
    file.exists(
      "/bin/bash"
    )
) {

  shell_command <- paste0(
    "source /etc/profile >/dev/null 2>&1 || true; ",
    "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; ",
    "module load imagemagick/7.1.0 >/dev/null 2>&1 || ",
    "module load imagemagick >/dev/null 2>&1 || ",
    "module load ImageMagick >/dev/null 2>&1 || true; ",
    "if command -v magick >/dev/null 2>&1; then ",
    "magick -density ",
    png_dpi,
    " '",
    pdf_file,
    "[0]' ",
    "-background white -alpha remove -alpha off '",
    png_file,
    "'; ",
    "elif command -v convert >/dev/null 2>&1; then ",
    "convert -density ",
    png_dpi,
    " '",
    pdf_file,
    "[0]' ",
    "-background white -alpha remove -alpha off '",
    png_file,
    "'; ",
    "else exit 127; fi"
  )

  shell_output <- tryCatch(
    system2(
      "/bin/bash",
      args = c(
        "-lc",
        shQuote(
          shell_command
        )
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      structure(
        conditionMessage(
          e
        ),
        status = 999L
      )
    }
  )

  shell_status <- attr(
    shell_output,
    "status"
  )

  if (
    is.null(
      shell_status
    )
  ) {
    shell_status <- 0L
  }

  Sys.sleep(
    1
  )

  if (
    shell_status ==
      0L &&
      png_is_valid(
        png_file
      )
  ) {

    png_method <- "ImageMagick via module/login shell"

    message(
      "PNG method succeeded: ",
      png_method
    )

  } else {

    message(
      "ImageMagick module fallback failed:\n",
      paste(
        shell_output,
        collapse = "\n"
      )
    )
  }
}

png_created <- !is.na(
  png_method
) &&
  png_is_valid(
    png_file
  )

if (
  !png_created
) {

  png_method <- NA_character_

  warning(
    "The complete 2011-2020 model and vector PDF were created, but this ",
    "Hopper session has no working headless PNG renderer. The PDF remains ",
    "the publication-quality master:\n",
    pdf_file
  )

  if (
    require_png
  ) {
    stop(
      "require_png = TRUE, but PNG creation failed."
    )
  }
}

# ------------------------------------------------------------------------------
# 16. Final verification — PNG is optional unless require_png = TRUE
# ------------------------------------------------------------------------------

required_outputs <- c(
  pdf_file,
  summary_qc_file,
  year_qc_file,
  parent_qc_file,
  plot_data_file,
  gasoline_output_file,
  diesel_output_file,
  preflight_inventory_file
)

if (
  require_png
) {
  required_outputs <- c(
    required_outputs,
    png_file
  )
}

missing_outputs <- required_outputs[
  !file.exists(
    required_outputs
  )
]

empty_outputs <- required_outputs[
  file.exists(
    required_outputs
  ) &
    (
      !is.finite(
        file.info(
          required_outputs
        )$size
      ) |
        file.info(
          required_outputs
        )$size <=
        0
    )
]

if (
  length(
    missing_outputs
  ) >
    0
) {
  stop(
    "Missing required final outputs:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}

if (
  length(
    empty_outputs
  ) >
    0
) {
  stop(
    "Empty required final outputs:\n",
    paste(
      empty_outputs,
      collapse = "\n"
    )
  )
}

completion_text <- c(
  paste0(
    "Script version: ",
    SCRIPT_VERSION
  ),
  paste0(
    "Completed years: ",
    paste(
      completed_years,
      collapse = ", "
    )
  ),
  "Sources: ONR gasoline; NRD diesel",
  paste0(
    "Wind method: eight-neighbor upwind angular weighting plus speed factor; U0=",
    U0_ms,
    ", beta=",
    beta_wind
  ),
  paste0(
    "Vector PDF: ",
    pdf_file
  ),
  paste0(
    "PNG created: ",
    png_created
  ),
  paste0(
    "PNG method: ",
    if (
      png_created
    ) {
      png_method
    } else {
      "not available in this Hopper session"
    }
  ),
  paste0(
    "Completed at: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z"
    )
  )
)

writeLines(
  completion_text,
  con = completion_file
)

all_outputs <- c(
  required_outputs,
  completion_file
)

if (
  png_created &&
    !png_file %in%
    all_outputs
) {
  all_outputs <- c(
    all_outputs,
    png_file
  )
}

output_summary <- data.table(
  file = all_outputs,
  size_mb = round(
    file.info(
      all_outputs
    )$size /
      1024^2,
    3
  )
)

cat(
  "\n============================================================\n",
  "FIGURE 4-4 DIRECTION+SPEED STRICT 2011-2020 RUN COMPLETED\n",
  "============================================================\n",
  "Completed years: ",
  paste(
    completed_years,
    collapse = ", "
  ),
  "\nVector PDF: ",
  pdf_file,
  "\nPNG created: ",
  png_created,
  "\nPNG method: ",
  if (
    png_created
  ) {
    png_method
  } else {
    "not available; use the vector PDF"
  },
  "\n",
  sep = ""
)

print(
  output_summary
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\nCompletion marker:\n",
  completion_file,
  "\n",
  sep = ""
)
