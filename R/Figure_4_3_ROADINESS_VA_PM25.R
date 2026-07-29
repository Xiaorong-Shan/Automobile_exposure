#!/usr/bin/env Rscript

# ==============================================================================
# Figures 4-2 and 4-3 -- final combined renderer
# Roadiness-based 1-km redistribution of gasoline- and diesel-attributable
# CMAQ-ISAM PM2.5 in Virginia, averaged over 2011-2020
#
# Confirmed source definitions:
#   ONR = gasoline vehicles
#   NRD = diesel vehicles
#
# IMPORTANT CORRECTION RELATIVE TO THE EARLIER SCRIPT:
#   - Gasoline roadiness uses highway + local-road LENGTH influence.
#   - Diesel roadiness uses highway-road LENGTH influence.
#   - Length-based roadiness is used because the manuscript defines roadiness
#     using road length and inverse-distance weighting, rather than segment count.
#
# Both figures use the same six-panel order:
#   (a) Gasoline: native 12-km CMAQ concentration
#   (b) Gasoline: roadiness-downscaled 1-km concentration
#   (c) Gasoline: 1-km minus native-parent concentration
#   (d) Diesel:   native 12-km CMAQ concentration
#   (e) Diesel:   roadiness-downscaled 1-km concentration
#   (f) Diesel:   1-km minus native-parent concentration
#
# Figure 4-2 shows the full Virginia domain.
# Figure 4-3 shows the Northern Virginia detail window.
#
# The native 12-km values are repeated over their assigned 1-km child cells only
# for visual comparison. Downscaling preserves the mean concentration among the
# modeled 1-km child cells assigned to each parent CMAQ grid cell.
#
# INPUTS:
#   /scratch/xshan2/R_Code/Automobiles/FIGURES/NATIVE_CMAQ_REBUILT/
#       DECADE_ONR_mean_2011_2020.rds
#       DECADE_NRD_mean_2011_2020.rds
#
#   /scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/VA/
#       roadiness_1km_hw_loc_VA.fst
#
# FIGURE OUTPUTS:
#   Figure_4_2_roadiness_VA_full_2011_2020.pdf
#   Figure_4_2_roadiness_VA_full_2011_2020.png
#   Figure_4_3_roadiness_NOVA_zoom_2011_2020.pdf
#   Figure_4_3_roadiness_NOVA_zoom_2011_2020.png
#
# SHARED DATA AND QC OUTPUTS:
#   Roadiness_QC_summary.csv
#   Roadiness_parent_QC.csv
#   Roadiness_plot_data.rds
#   VA_1km_gasoline_roadiness_2011_2020_mean.fst
#   VA_1km_diesel_roadiness_2011_2020_mean.fst
#
# HOPPER:
#   The script always creates a vector PDF with base R and then attempts several
#   independent headless PNG methods. It reports success only after every output
#   actually exists.
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
  library(viridis)
})

# ------------------------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------------------------

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

gasoline_file <- file.path(
  cmaq_dir,
  "DECADE_ONR_mean_2011_2020.rds"
)

diesel_file <- file.path(
  cmaq_dir,
  "DECADE_NRD_mean_2011_2020.rds"
)

roadiness_file <- paste0(
  "/scratch/xshan2/R_Code/disperseR/Auto/",
  "roadiness_2017/VA/roadiness_1km_hw_loc_VA.fst"
)

out_dir <- file.path(
  cmaq_dir,
  "FIGURES_4_2_4_3_ROADINESS_VA_AND_NOVA"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

statewide_pdf_file <- file.path(
  out_dir,
  "Figure_4_2_roadiness_VA_full_2011_2020.pdf"
)

statewide_png_file <- file.path(
  out_dir,
  "Figure_4_2_roadiness_VA_full_2011_2020.png"
)

nova_pdf_file <- file.path(
  out_dir,
  "Figure_4_3_roadiness_NOVA_zoom_2011_2020.pdf"
)

nova_png_file <- file.path(
  out_dir,
  "Figure_4_3_roadiness_NOVA_zoom_2011_2020.png"
)

# The legacy renderer helpers below use pdf_file and png_file as the active
# output targets. Start with Figure 4-3, then switch them to Figure 4-2 later.
pdf_file <- nova_pdf_file
png_file <- nova_png_file

summary_qc_file <- file.path(
  out_dir,
  "Roadiness_QC_summary.csv"
)

parent_qc_file <- file.path(
  out_dir,
  "Roadiness_parent_QC.csv"
)

plot_data_file <- file.path(
  out_dir,
  "Roadiness_plot_data.rds"
)

gasoline_output_file <- file.path(
  out_dir,
  "VA_1km_gasoline_roadiness_2011_2020_mean.fst"
)

diesel_output_file <- file.path(
  out_dir,
  "VA_1km_diesel_roadiness_2011_2020_mean.fst"
)

required_inputs <- c(
  gasoline_file,
  diesel_file,
  roadiness_file
)

missing_inputs <- required_inputs[
  !file.exists(required_inputs)
]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required input files:\n",
    paste(
      missing_inputs,
      collapse = "\n"
    )
  )
}

# Display settings
# Fixed absolute PM2.5 legend range for all native and downscaled panels.
# Values above 1.0 ug/m3 are displayed at the upper color, but the data are unchanged.
absolute_upper_limit_fixed <- 1.0
difference_upper_quantile <- 0.995

# Separate dimensions let the full-state and zoomed figures fill their panels.
# PDFs remain vector originals; 600-dpi PNGs are intended for insertion in Word.
statewide_width_in <- 14.4
statewide_height_in <- 8.0
nova_width_in <- 14.4
nova_height_in <- 8.8
png_dpi <- 600

# Active dimensions for the first-rendered Northern Virginia figure.
fig_width_in <- nova_width_in
fig_height_in <- nova_height_in

# Numerical tolerance used only for QC
small_value <- 1e-12

# ------------------------------------------------------------------------------
# 2. Read the two native CMAQ decade-mean surfaces
# ------------------------------------------------------------------------------

read_cmaq_surface <- function(
    file_path,
    source_code,
    source_name) {

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
    names(dt)
  )

  if (length(missing_columns) > 0) {
    stop(
      "CMAQ file is missing required columns:\n",
      file_path,
      "\nMissing: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  if ("n_years" %in% names(dt)) {

    out <- dt[, .(
      lon = as.numeric(x),
      lat = as.numeric(y),
      concentration = as.numeric(value),
      n_years = as.integer(n_years)
    )]

  } else {

    out <- dt[, .(
      lon = as.numeric(x),
      lat = as.numeric(y),
      concentration = as.numeric(value),
      n_years = NA_integer_
    )]
  }

  out <- out[
    is.finite(lon) &
      is.finite(lat) &
      is.finite(concentration)
  ]

  if (nrow(out) == 0) {
    stop(
      "No finite CMAQ values were found in:\n",
      file_path
    )
  }

  if (
    anyDuplicated(
      out[, .(
        lon,
        lat
      )]
    ) > 0
  ) {
    stop(
      "Duplicate CMAQ coordinates were found in:\n",
      file_path
    )
  }

  out[, source_code := source_code]
  out[, source_name := source_name]

  out
}

gasoline_cmaq <- read_cmaq_surface(
  file_path = gasoline_file,
  source_code = "ONR",
  source_name = "Gasoline vehicles"
)

diesel_cmaq <- read_cmaq_surface(
  file_path = diesel_file,
  source_code = "NRD",
  source_name = "Diesel vehicles"
)

gasoline_coordinates <- gasoline_cmaq[
  order(
    lon,
    lat
  ),
  .(
    lon,
    lat
  )
]

diesel_coordinates <- diesel_cmaq[
  order(
    lon,
    lat
  ),
  .(
    lon,
    lat
  )
]

if (!identical(
  gasoline_coordinates,
  diesel_coordinates
)) {
  stop(
    "Gasoline and diesel CMAQ surfaces do not use identical coordinates."
  )
}

# ------------------------------------------------------------------------------
# 3. Read roadiness and construct source-specific length-based metrics
# ------------------------------------------------------------------------------

road_columns <- c(
  "x",
  "y",
  "lengm.dist2_hw",
  "lengm.dist2_loc"
)

road_dt <- read_fst(
  roadiness_file,
  columns = road_columns,
  as.data.table = TRUE
)

missing_road_columns <- setdiff(
  road_columns,
  names(road_dt)
)

if (length(missing_road_columns) > 0) {
  stop(
    "Roadiness FST is missing required columns:\n",
    paste(
      missing_road_columns,
      collapse = ", "
    )
  )
}

road_dt <- road_dt[, .(
  x = as.numeric(x),
  y = as.numeric(y),
  road_highway = as.numeric(lengm.dist2_hw),
  road_local = as.numeric(lengm.dist2_loc)
)]

road_dt[
  !is.finite(road_highway) |
    is.na(road_highway) |
    road_highway < 0,
  road_highway := 0
]

road_dt[
  !is.finite(road_local) |
    is.na(road_local) |
    road_local < 0,
  road_local := 0
]

if (
  anyDuplicated(
    road_dt[, .(
      x,
      y
    )]
  ) > 0
) {
  stop(
    "Duplicate 1-km x-y coordinates were found in the roadiness FST."
  )
}

# Roadiness-only source assumptions:
# - Gasoline vehicles use both highway and local-road length influence.
# - Diesel vehicles use highway-road length influence.
road_dt[, road_gasoline :=
          road_highway + road_local]

road_dt[, road_diesel :=
          road_highway]

# ------------------------------------------------------------------------------
# 4. Virginia boundary and clipping of the 1-km template
# ------------------------------------------------------------------------------

states_ll <- USAboundaries::us_states(
  resolution = "low"
)

virginia_ll <- states_ll[
  states_ll$state_abbr == "VA",
]

if (nrow(virginia_ll) != 1) {
  stop(
    "Virginia was not uniquely identified in USAboundaries."
  )
}

virginia_lcc <- st_transform(
  virginia_ll,
  crs = st_crs(p4s)
)

road_points_lcc <- st_as_sf(
  road_dt,
  coords = c(
    "x",
    "y"
  ),
  crs = st_crs(p4s),
  remove = FALSE
)

inside_va <- lengths(
  st_intersects(
    road_points_lcc,
    virginia_lcc
  )
) > 0

road_points_lcc <- road_points_lcc[
  inside_va,
]

road_dt <- as.data.table(
  st_drop_geometry(
    road_points_lcc
  )
)

if (nrow(road_dt) == 0) {
  stop(
    "No 1-km roadiness grid-cell centers were retained inside Virginia."
  )
}

cat(
  "\nVirginia 1-km cells retained: ",
  format(
    nrow(road_dt),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 5. Build the nearby native CMAQ grid and map 1-km cells to parent centers
# ------------------------------------------------------------------------------

va_bbox_ll <- st_bbox(
  virginia_ll
)

# Add a generous geographic margin so every Virginia 1-km cell has nearby
# native CMAQ centers available for nearest-center assignment.
cmaq_margin_deg <- 1.5

cmaq_grid <- unique(
  gasoline_cmaq[
    lon >= as.numeric(va_bbox_ll["xmin"]) - cmaq_margin_deg &
      lon <= as.numeric(va_bbox_ll["xmax"]) + cmaq_margin_deg &
      lat >= as.numeric(va_bbox_ll["ymin"]) - cmaq_margin_deg &
      lat <= as.numeric(va_bbox_ll["ymax"]) + cmaq_margin_deg,
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

if (nrow(cmaq_grid) == 0) {
  stop(
    "No native CMAQ grid centers were retained around Virginia."
  )
}

cmaq_points_ll <- st_as_sf(
  cmaq_grid,
  coords = c(
    "lon",
    "lat"
  ),
  crs = 4326,
  remove = FALSE
)

cmaq_points_lcc <- st_transform(
  cmaq_points_ll,
  crs = st_crs(p4s)
)

parent_index <- st_nearest_feature(
  road_points_lcc,
  cmaq_points_lcc
)

if (
  length(parent_index) !=
    nrow(road_dt)
) {
  stop(
    "Nearest-parent index length does not match the 1-km roadiness grid."
  )
}

nearest_parent_distance_m <- as.numeric(
  st_distance(
    road_points_lcc,
    cmaq_points_lcc[parent_index, ],
    by_element = TRUE
  )
)

if (
  any(!is.finite(nearest_parent_distance_m))
) {
  stop(
    "Non-finite nearest-parent distances were produced."
  )
}

road_dt[, parent_id :=
          cmaq_grid$parent_id[parent_index]]

road_dt[, parent_distance_m :=
          nearest_parent_distance_m]

cat(
  "Nearest CMAQ-center distance, median: ",
  round(
    median(nearest_parent_distance_m),
    1
  ),
  " m; maximum: ",
  round(
    max(nearest_parent_distance_m),
    1
  ),
  " m\n",
  sep = ""
)

# A regular 12-km grid should normally place all child centers within roughly
# 8.5 km of the nearest parent center. This is a warning, not a forced stop,
# because edge geometry and projection can increase the maximum slightly.
if (
  max(
    nearest_parent_distance_m,
    na.rm = TRUE
  ) > 20000
) {
  warning(
    "Some 1-km cells are more than 20 km from the assigned CMAQ center. ",
    "Inspect the CRS and grid alignment before using the result."
  )
}

# ------------------------------------------------------------------------------
# 6. Attach native source concentrations to the parent grid
# ------------------------------------------------------------------------------

attach_parent_values <- function(
    source_dt,
    source_code,
    source_name) {

  parent_values <- merge(
    cmaq_grid,
    source_dt[, .(
      lon,
      lat,
      pm12 = concentration
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
      parent_values$pm12
    )
  ) {
    stop(
      "Some nearby parent CMAQ cells are missing ",
      source_name,
      " concentrations."
    )
  }

  parent_values[, source_code := source_code]
  parent_values[, source_name := source_name]

  parent_values
}

gasoline_parent_values <- attach_parent_values(
  source_dt = gasoline_cmaq,
  source_code = "ONR",
  source_name = "Gasoline vehicles"
)

diesel_parent_values <- attach_parent_values(
  source_dt = diesel_cmaq,
  source_code = "NRD",
  source_name = "Diesel vehicles"
)

# ------------------------------------------------------------------------------
# 7. Mean-preserving roadiness redistribution
# ------------------------------------------------------------------------------

downscale_one_source <- function(
    parent_values,
    roadiness_column,
    source_code,
    source_name) {

  dt <- copy(
    road_dt[, .(
      x,
      y,
      parent_id,
      parent_distance_m
    )]
  )

  dt[, pm12 :=
       parent_values$pm12[parent_id]]

  dt[, roadiness :=
       as.numeric(
         road_dt[[roadiness_column]]
       )]

  dt[
    !is.finite(roadiness) |
      is.na(roadiness) |
      roadiness < 0,
    roadiness := 0
  ]

  # Normalize roadiness within each assigned parent grid cell.
  dt[, parent_roadiness_mean :=
       mean(
         roadiness,
         na.rm = TRUE
       ),
     by = parent_id]

  dt[, weight :=
       fifelse(
         is.finite(parent_roadiness_mean) &
           parent_roadiness_mean > 0,
         roadiness / parent_roadiness_mean,
         1
       )]

  if (
    any(!is.finite(dt$weight))
  ) {
    stop(
      "Non-finite downscaling weights were produced for ",
      source_name,
      "."
    )
  }

  dt[, pm1 :=
       pm12 * weight]

  dt[, delta :=
       pm1 - pm12]

  dt[, source_code := source_code]
  dt[, source_name := source_name]

  dt
}

gasoline_1km <- downscale_one_source(
  parent_values = gasoline_parent_values,
  roadiness_column = "road_gasoline",
  source_code = "ONR",
  source_name = "Gasoline vehicles"
)

diesel_1km <- downscale_one_source(
  parent_values = diesel_parent_values,
  roadiness_column = "road_diesel",
  source_code = "NRD",
  source_name = "Diesel vehicles"
)

plot_data <- rbindlist(
  list(
    gasoline_1km,
    diesel_1km
  ),
  use.names = TRUE
)

plot_data[, source_name :=
            factor(
              source_name,
              levels = c(
                "Gasoline vehicles",
                "Diesel vehicles"
              )
            )]

# Save the model-ready source-specific 1-km products.
write_fst(
  gasoline_1km,
  gasoline_output_file
)

write_fst(
  diesel_1km,
  diesel_output_file
)

# ------------------------------------------------------------------------------
# 8. Quantitative QC and mass-conservation diagnostics
# ------------------------------------------------------------------------------

parent_qc <- plot_data[, {

  parent_pm12 <- first(pm12)
  parent_pm1_mean <- mean(
    pm1,
    na.rm = TRUE
  )

  parent_cv <- if (
    .N > 1 &&
      is.finite(parent_pm1_mean) &&
      parent_pm1_mean > 0
  ) {
    sd(
      pm1,
      na.rm = TRUE
    ) / parent_pm1_mean
  } else {
    0
  }

  list(
    n_1km_cells = .N,
    pm12 = parent_pm12,
    mean_pm1 = parent_pm1_mean,
    absolute_mass_error = parent_pm1_mean - parent_pm12,
    relative_mass_error = if (
      is.finite(parent_pm12) &&
        abs(parent_pm12) > small_value
    ) {
      (
        parent_pm1_mean -
          parent_pm12
      ) / parent_pm12
    } else {
      NA_real_
    },
    within_parent_cv = parent_cv,
    within_parent_p05 = as.numeric(
      quantile(
        pm1,
        0.05,
        na.rm = TRUE
      )
    ),
    within_parent_p95 = as.numeric(
      quantile(
        pm1,
        0.95,
        na.rm = TRUE
      )
    )
  )

}, by = .(
  source_code,
  source_name,
  parent_id
)]

summary_qc <- plot_data[, .(
  n_1km_cells = .N,
  n_parent_cells = uniqueN(parent_id),

  mean_native_child_representation = mean(
    pm12,
    na.rm = TRUE
  ),

  mean_downscaled = mean(
    pm1,
    na.rm = TRUE
  ),

  median_downscaled = median(
    pm1,
    na.rm = TRUE
  ),

  downscaled_p05 = as.numeric(
    quantile(
      pm1,
      0.05,
      na.rm = TRUE
    )
  ),

  downscaled_p95 = as.numeric(
    quantile(
      pm1,
      0.95,
      na.rm = TRUE
    )
  ),

  downscaled_p99 = as.numeric(
    quantile(
      pm1,
      0.99,
      na.rm = TRUE
    )
  ),

  delta_p05 = as.numeric(
    quantile(
      delta,
      0.05,
      na.rm = TRUE
    )
  ),

  delta_p95 = as.numeric(
    quantile(
      delta,
      0.95,
      na.rm = TRUE
    )
  ),

  maximum_positive_delta = max(
    delta,
    na.rm = TRUE
  ),

  minimum_negative_delta = min(
    delta,
    na.rm = TRUE
  ),

  percent_at_least_50_percent_above_parent =
    100 * mean(
      pm1 >=
        1.5 * pm12,
      na.rm = TRUE
    ),

  percent_at_least_50_percent_below_parent =
    100 * mean(
      pm1 <=
        0.5 * pm12,
      na.rm = TRUE
    ),

  median_roadiness_weight = median(
    weight,
    na.rm = TRUE
  ),

  roadiness_weight_p95 = as.numeric(
    quantile(
      weight,
      0.95,
      na.rm = TRUE
    )
  ),

  roadiness_weight_p99 = as.numeric(
    quantile(
      weight,
      0.99,
      na.rm = TRUE
    )
  )

), by = .(
  source_code,
  source_name
)]

parent_summary <- parent_qc[, .(
  median_within_parent_cv = median(
    within_parent_cv,
    na.rm = TRUE
  ),

  p95_within_parent_cv = as.numeric(
    quantile(
      within_parent_cv,
      0.95,
      na.rm = TRUE
    )
  ),

  maximum_absolute_parent_mean_error = max(
    abs(
      absolute_mass_error
    ),
    na.rm = TRUE
  ),

  maximum_absolute_relative_parent_mean_error = max(
    abs(
      relative_mass_error
    ),
    na.rm = TRUE
  )

), by = .(
  source_code,
  source_name
)]

summary_qc <- merge(
  summary_qc,
  parent_summary,
  by = c(
    "source_code",
    "source_name"
  ),
  all.x = TRUE,
  sort = FALSE
)

summary_qc[, median_parent_assignment_distance_m :=
             median(
               plot_data[
                 source_code ==
                   .BY$source_code,
                 parent_distance_m
               ],
               na.rm = TRUE
             ),
           by = source_code]

summary_qc[, maximum_parent_assignment_distance_m :=
             max(
               plot_data[
                 source_code ==
                   .BY$source_code,
                 parent_distance_m
               ],
               na.rm = TRUE
             ),
           by = source_code]

fwrite(
  summary_qc,
  summary_qc_file
)

fwrite(
  parent_qc,
  parent_qc_file
)

saveRDS(
  list(
    plot_data = plot_data,
    summary_qc = summary_qc,
    parent_qc = parent_qc,
    gasoline_roadiness_definition =
      "lengm.dist2_hw + lengm.dist2_loc",
    diesel_roadiness_definition =
      "lengm.dist2_hw"
  ),
  plot_data_file
)

cat(
  "\n========================================\n",
  "SHARED FIGURES 4-2 AND 4-3 QC SUMMARY\n",
  "========================================\n"
)

print(
  summary_qc
)

# ------------------------------------------------------------------------------
# 9. Shared map scales
# ------------------------------------------------------------------------------

absolute_values <- c(
  plot_data$pm12,
  plot_data$pm1
)

# Fixed common absolute-concentration scale used in both figures.
absolute_upper_limit <- absolute_upper_limit_fixed
absolute_breaks <- seq(
  0,
  absolute_upper_limit,
  by = 0.2
)

# Diagnostic only: report the fraction of plotted values that will be
# displayed at the upper color because they exceed 1.0 ug/m3.
percent_absolute_values_above_legend <- 100 * mean(
  absolute_values > absolute_upper_limit,
  na.rm = TRUE
)

cat(
  "Percent of native/downscaled PM2.5 values above the fixed 1.0 ug/m3 display limit: ",
  round(percent_absolute_values_above_legend, 3),
  "%\n",
  sep = ""
)

difference_upper_limit <- as.numeric(
  quantile(
    abs(
      plot_data$delta
    ),
    probs = difference_upper_quantile,
    na.rm = TRUE
  )
)

if (
  !is.finite(difference_upper_limit) ||
    difference_upper_limit <= 0
) {
  difference_upper_limit <- max(
    abs(
      plot_data$delta
    ),
    na.rm = TRUE
  )
}

difference_breaks <- pretty(
  c(
    -difference_upper_limit,
    difference_upper_limit
  ),
  n = 5
)

difference_breaks <- difference_breaks[
  difference_breaks >=
    -difference_upper_limit &
    difference_breaks <=
    difference_upper_limit
]

if (
  length(
    difference_breaks
  ) < 3
) {
  difference_breaks <- seq(
    -difference_upper_limit,
    difference_upper_limit,
    length.out = 5
  )
}

cat(
  "\nAbsolute display upper limit: ",
  signif(
    absolute_upper_limit,
    5
  ),
  " ug/m3\n",
  "Difference display limits: +/-",
  signif(
    difference_upper_limit,
    5
  ),
  " ug/m3\n",
  sep = ""
)

# Common Northern Virginia zoom window for all six panels. Keeping the extent
# in longitude/latitude makes the selected section transparent and easy to
# change if a different Virginia subregion is requested later.
#
# Current window: Northern Virginia, including the Fairfax-Arlington-Prince
# William-Loudoun corridor and its immediate surroundings.
zoom_bbox_ll <- st_bbox(
  c(
    xmin = -78.25,
    ymin = 38.25,
    xmax = -76.80,
    ymax = 39.50
  ),
  crs = st_crs(
    4326
  )
)

zoom_bbox_lcc <- st_bbox(
  st_transform(
    st_as_sfc(
      zoom_bbox_ll
    ),
    st_crs(
      virginia_lcc
    )
  )
)

zoom_padding_m <- 2500

nova_map_xlim <- c(
  as.numeric(
    zoom_bbox_lcc["xmin"]
  ) - zoom_padding_m,
  as.numeric(
    zoom_bbox_lcc["xmax"]
  ) + zoom_padding_m
)

nova_map_ylim <- c(
  as.numeric(
    zoom_bbox_lcc["ymin"]
  ) - zoom_padding_m,
  as.numeric(
    zoom_bbox_lcc["ymax"]
  ) + zoom_padding_m
)

# Active map bounds used by the existing panel-building calls.
map_xlim <- nova_map_xlim
map_ylim <- nova_map_ylim

cat(
  "\nFigure zoom window (longitude/latitude):\n",
  "  xmin = ",
  zoom_bbox_ll["xmin"],
  "\n  ymin = ",
  zoom_bbox_ll["ymin"],
  "\n  xmax = ",
  zoom_bbox_ll["xmax"],
  "\n  ymax = ",
  zoom_bbox_ll["ymax"],
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 10. Publication-style map panels
# ------------------------------------------------------------------------------

make_absolute_panel <- function(
    data,
    value_column,
    map_xlim_use = map_xlim,
    map_ylim_use = map_ylim,
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
      linewidth = 0.52
    ) +

    coord_sf(
      crs = st_crs(p4s),
      xlim = map_xlim_use,
      ylim = map_ylim_use,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +

    scale_fill_viridis_c(
      option = "E",
      direction = 1,
      limits = c(
        0,
        absolute_upper_limit
      ),
      breaks = absolute_breaks,
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
          92,
          "mm"
        ),
        barheight = grid::unit(
          4.6,
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
        linewidth = 0.42
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
        size = 12.0
      ),

      legend.text = element_text(
        size = 10.5
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
        2,
        3,
        3,
        3,
        unit = "pt"
      )
    )
}

make_difference_panel <- function(
    data,
    map_xlim_use = map_xlim,
    map_ylim_use = map_ylim,
    show_legend = FALSE) {

  ggplot() +

    geom_raster(
      data = data,
      aes(
        x = x,
        y = y,
        fill = delta
      ),
      interpolate = FALSE
    ) +

    geom_sf(
      data = virginia_lcc,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey25",
      linewidth = 0.52
    ) +

    coord_sf(
      crs = st_crs(p4s),
      xlim = map_xlim_use,
      ylim = map_ylim_use,
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
        -difference_upper_limit,
        difference_upper_limit
      ),
      breaks = difference_breaks,
      labels = scales::label_number(
        accuracy = 0.1
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
          56,
          "mm"
        ),
        barheight = grid::unit(
          4.6,
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
        linewidth = 0.42
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
        size = 12.0
      ),

      legend.text = element_text(
        size = 10.5
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
        2,
        3,
        3,
        3,
        unit = "pt"
      )
    )
}

gasoline_plot_data <- plot_data[
  source_code == "ONR"
]

diesel_plot_data <- plot_data[
  source_code == "NRD"
]

p_a <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm12",
  show_legend = FALSE
)

p_b <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm1",
  show_legend = FALSE
)

p_c <- make_difference_panel(
  data = gasoline_plot_data,
  show_legend = FALSE
)

p_d <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm12",
  show_legend = FALSE
)

p_e <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm1",
  show_legend = FALSE
)

p_f <- make_difference_panel(
  data = diesel_plot_data,
  show_legend = FALSE
)

# Legend-only plots
p_absolute_legend <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm1",
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
    ) == 0
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
# 11. Composite figure drawing
# ------------------------------------------------------------------------------

draw_figure <- function() {

  grid::grid.newpage()

  outer_layout <- grid::grid.layout(
    nrow = 4,
    ncol = 5,
    widths = grid::unit.c(
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        0.14,
        "in"
      ),
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        0.14,
        "in"
      ),
      grid::unit(
        1,
        "null"
      )
    ),
    heights = grid::unit.c(
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        0.13,
        "in"
      ),
      grid::unit(
        1,
        "null"
      ),
      grid::unit(
        0.84,
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
          0.27,
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
        2.0,
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
        fontsize = 15.0
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
    3
  )

  draw_panel_with_tag(
    g_c,
    "(c)",
    1,
    5
  )

  draw_panel_with_tag(
    g_d,
    "(d)",
    3,
    1
  )

  draw_panel_with_tag(
    g_e,
    "(e)",
    3,
    3
  )

  draw_panel_with_tag(
    g_f,
    "(f)",
    3,
    5
  )

  # Shared absolute concentration legend under columns 1-2
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 4,
      layout.pos.col = 1:3
    )
  )

  grid::grid.draw(
    absolute_legend_grob
  )

  grid::popViewport()

  # Shared difference legend under column 3
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 4,
      layout.pos.col = 5
    )
  )

  grid::grid.draw(
    difference_legend_grob
  )

  grid::popViewport()

  grid::popViewport()
}

# ------------------------------------------------------------------------------
# 12. Save vector PDF
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
    )$size <= 1000
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
# 12A. Save the full-Virginia PDF before attempting either PNG
# ------------------------------------------------------------------------------
#
# The earlier combined script attempted the Northern Virginia PNG before it
# reached the full-state PDF section. A missing headless PNG renderer therefore
# stopped the script with only the NOVA PDF on disk. The full-state vector PDF
# is now created and verified here, before any PNG work begins.

nova_map_xlim_saved <- map_xlim
nova_map_ylim_saved <- map_ylim

va_bbox_lcc_early <- st_bbox(
  virginia_lcc
)

statewide_padding_m_early <- 10000

map_xlim <- c(
  as.numeric(
    va_bbox_lcc_early["xmin"]
  ) - statewide_padding_m_early,
  as.numeric(
    va_bbox_lcc_early["xmax"]
  ) + statewide_padding_m_early
)

map_ylim <- c(
  as.numeric(
    va_bbox_lcc_early["ymin"]
  ) - statewide_padding_m_early,
  as.numeric(
    va_bbox_lcc_early["ymax"]
  ) + statewide_padding_m_early
)

rebuild_active_figure_grobs <- function() {

  p_a <<- make_absolute_panel(
    data = gasoline_plot_data,
    value_column = "pm12",
    show_legend = FALSE
  )

  p_b <<- make_absolute_panel(
    data = gasoline_plot_data,
    value_column = "pm1",
    show_legend = FALSE
  )

  p_c <<- make_difference_panel(
    data = gasoline_plot_data,
    show_legend = FALSE
  )

  p_d <<- make_absolute_panel(
    data = diesel_plot_data,
    value_column = "pm12",
    show_legend = FALSE
  )

  p_e <<- make_absolute_panel(
    data = diesel_plot_data,
    value_column = "pm1",
    show_legend = FALSE
  )

  p_f <<- make_difference_panel(
    data = diesel_plot_data,
    show_legend = FALSE
  )

  p_absolute_legend <<- make_absolute_panel(
    data = gasoline_plot_data,
    value_column = "pm1",
    show_legend = TRUE
  ) +
    theme(
      panel.border = element_blank()
    )

  p_difference_legend <<- make_difference_panel(
    data = gasoline_plot_data,
    show_legend = TRUE
  ) +
    theme(
      panel.border = element_blank()
    )

  absolute_legend_grob <<- get_legend(
    p_absolute_legend
  )

  difference_legend_grob <<- get_legend(
    p_difference_legend
  )

  g_a <<- ggplotGrob(
    p_a
  )

  g_b <<- ggplotGrob(
    p_b
  )

  g_c <<- ggplotGrob(
    p_c
  )

  g_d <<- ggplotGrob(
    p_d
  )

  g_e <<- ggplotGrob(
    p_e
  )

  g_f <<- ggplotGrob(
    p_f
  )

  invisible(
    NULL
  )
}

rebuild_active_figure_grobs()

unlink(
  statewide_pdf_file,
  force = TRUE
)

grDevices::pdf(
  file = statewide_pdf_file,
  width = statewide_width_in,
  height = statewide_height_in,
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
    statewide_pdf_file
  ) ||
    !is.finite(
      file.info(
        statewide_pdf_file
      )$size
    ) ||
    file.info(
      statewide_pdf_file
    )$size <= 1000
) {
  stop(
    "Figure 4-2 full-Virginia PDF creation failed:\n",
    statewide_pdf_file
  )
}

cat(
  "\nFigure 4-2 full-Virginia vector PDF created before PNG rendering:\n",
  statewide_pdf_file,
  "\nSize: ",
  round(
    file.info(
      statewide_pdf_file
    )$size /
      1024^2,
    3
  ),
  " MB\n",
  sep = ""
)

# Restore and rebuild the NOVA grobs so the next section creates the correct
# Figure 4-3 PNG rather than accidentally rasterizing the statewide map.
map_xlim <- nova_map_xlim_saved
map_ylim <- nova_map_ylim_saved

rebuild_active_figure_grobs()

# ------------------------------------------------------------------------------
# 13. Hopper-safe PNG creation with independent fallbacks
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
    )$size > 1000
}

close_extra_devices <- function(
    start_device) {

  while (
    grDevices::dev.cur() > 1 &&
      grDevices::dev.cur() != start_device
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

# Method 1: base-R cairo PNG, when compiled into this R installation.
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

# Method 2: ragg, when installed.
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

# Method 3: Cairo package, only when installed.
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

# Method 4: base bitmap / Ghostscript.
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
    status != 0L
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
    status == 0L &&
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

# Method 5: Ghostscript PDF -> PNG.
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
    ) == 1 &&
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

# Method 6: pdftocairo.
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
    ) == 1 &&
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

# Method 7: pdftoppm.
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
    ) == 1 &&
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

# Method 8: ImageMagick on PATH.
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
    ) == 1 &&
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
    ) == 1 &&
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

# Final fallback: ImageMagick through the Hopper module/login environment.
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
    shell_status == 0L &&
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

if (
  is.na(
    png_method
  ) ||
    !png_is_valid(
      png_file
    )
) {
  stop(
    "The vector PDF was created successfully, but no available headless PNG ",
    "renderer/converter was found in this Hopper session.\n\n",
    "PDF:\n",
    pdf_file,
    "\n\n",
    "The script tried base cairo, ragg, Cairo, bitmap/Ghostscript, gs, ",
    "pdftocairo, pdftoppm, and ImageMagick."
  )
}

# ------------------------------------------------------------------------------
# 14. Final verification
# ------------------------------------------------------------------------------

expected_outputs <- c(
  pdf_file,
  png_file,
  statewide_pdf_file,
  summary_qc_file,
  parent_qc_file,
  plot_data_file,
  gasoline_output_file,
  diesel_output_file
)

missing_outputs <- expected_outputs[
  !file.exists(
    expected_outputs
  )
]

empty_outputs <- expected_outputs[
  file.exists(
    expected_outputs
  ) &
    (
      !is.finite(
        file.info(
          expected_outputs
        )$size
      ) |
        file.info(
          expected_outputs
        )$size <= 0
    )
]

if (
  length(
    missing_outputs
  ) > 0
) {
  stop(
    "Missing expected outputs:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}

if (
  length(
    empty_outputs
  ) > 0
) {
  stop(
    "Empty expected outputs:\n",
    paste(
      empty_outputs,
      collapse = "\n"
    )
  )
}

output_summary <- data.table(
  file = expected_outputs,
  size_mb = round(
    file.info(
      expected_outputs
    )$size /
      1024^2,
    3
  )
)

cat(
  "\n========================================\n",
  "FIGURE 4-3 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "PNG method used: ",
  png_method,
  "\n",
  sep = ""
)

print(
  output_summary
)

cat(
  "\nOutput directory:\n",
  out_dir,
  "\n",
  sep = ""
)

# ==============================================================================
# 15. Figure 4-2: full Virginia domain
# ==============================================================================
#
# The calculation, scales, panel order, and legend breaks are identical to
# Figure 4-3. Only the map extent and output canvas height change. Reusing the
# already validated plot_data prevents the statewide and zoom figures from
# diverging because of separate calculations.

nova_pdf_file <- pdf_file
nova_png_file <- png_file
nova_png_method <- png_method

statewide_pdf_file <- file.path(
  out_dir,
  "Figure_4_2_roadiness_VA_full_2011_2020.pdf"
)

statewide_png_file <- file.path(
  out_dir,
  "Figure_4_2_roadiness_VA_full_2011_2020.png"
)

va_bbox_lcc <- st_bbox(
  virginia_lcc
)

statewide_padding_m <- 10000

map_xlim <- c(
  as.numeric(
    va_bbox_lcc["xmin"]
  ) - statewide_padding_m,
  as.numeric(
    va_bbox_lcc["xmax"]
  ) + statewide_padding_m
)

map_ylim <- c(
  as.numeric(
    va_bbox_lcc["ymin"]
  ) - statewide_padding_m,
  as.numeric(
    va_bbox_lcc["ymax"]
  ) + statewide_padding_m
)

cat(
  "\nRebuilding the same six panels for the full Virginia domain.\n"
)

# Rebuild all six panels after changing map_xlim and map_ylim. The plotting
# functions evaluate these bounds when each ggplot object is created.
p_a <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm12",
  show_legend = FALSE
)

p_b <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm1",
  show_legend = FALSE
)

p_c <- make_difference_panel(
  data = gasoline_plot_data,
  show_legend = FALSE
)

p_d <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm12",
  show_legend = FALSE
)

p_e <- make_absolute_panel(
  data = diesel_plot_data,
  value_column = "pm1",
  show_legend = FALSE
)

p_f <- make_difference_panel(
  data = diesel_plot_data,
  show_legend = FALSE
)

p_absolute_legend <- make_absolute_panel(
  data = gasoline_plot_data,
  value_column = "pm1",
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

# The statewide maps are wider than the Northern Virginia window. A slightly
# shorter canvas lets the Virginia boundary fill the panels without creating
# unnecessary vertical white space.
pdf_file <- statewide_pdf_file
png_file <- statewide_png_file
fig_width_in <- 14.4
fig_height_in <- 8.0

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
    )$size <= 1000
) {
  stop(
    "Figure 4-2 vector PDF creation failed:\n",
    pdf_file
  )
}

cat(
  "\nFigure 4-2 vector PDF created:\n",
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

# Reuse the tested Hopper-safe device helpers from Figure 4-3.
statewide_png_method <- NA_character_

if (
  isTRUE(
    capabilities(
      "cairo"
    )
  )
) {

  ok <- try_png_device(
    "Figure 4-2 base grDevices::png(type='cairo-png')",
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
    statewide_png_method <- "base cairo-png"
  }
}

if (
  is.na(
    statewide_png_method
  ) &&
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
) {

  ok <- try_png_device(
    "Figure 4-2 ragg::agg_png",
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
    statewide_png_method <- "ragg"
  }
}

if (
  is.na(
    statewide_png_method
  ) &&
    requireNamespace(
      "Cairo",
      quietly = TRUE
    )
) {

  ok <- try_png_device(
    "Figure 4-2 Cairo::CairoPNG",
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
    statewide_png_method <- "Cairo package"
  }
}

if (
  is.na(
    statewide_png_method
  )
) {

  ok <- try_png_device(
    "Figure 4-2 grDevices::bitmap(type='png16m')",
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
    statewide_png_method <- "base bitmap / Ghostscript"
  }
}

# Conversion from the vector PDF is retained as a fallback when no direct
# headless PNG device is available in the Hopper R installation.
if (
  is.na(
    statewide_png_method
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
    ) == 1 &&
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
      method_name = "Figure 4-2 pdftocairo"
    )

    if (
      ok
    ) {
      statewide_png_method <- "pdftocairo"
    }
  }
}

if (
  is.na(
    statewide_png_method
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
    ) == 1 &&
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
      method_name = "Figure 4-2 Ghostscript"
    )

    if (
      ok
    ) {
      statewide_png_method <- "Ghostscript"
    }
  }
}

if (
  is.na(
    statewide_png_method
  ) ||
    !png_is_valid(
      png_file
    )
) {
  stop(
    "Figure 4-2 vector PDF was created, but the 600-dpi PNG could not be ",
    "created by any available Hopper renderer.\nPDF:\n",
    pdf_file
  )
}

combined_figure_outputs <- c(
  statewide_pdf_file,
  statewide_png_file,
  nova_pdf_file,
  nova_png_file
)

missing_figure_outputs <- combined_figure_outputs[
  !file.exists(
    combined_figure_outputs
  ) |
    !is.finite(
      file.info(
        combined_figure_outputs
      )$size
    ) |
    file.info(
      combined_figure_outputs
    )$size <= 1000
]

if (
  length(
    missing_figure_outputs
  ) > 0
) {
  stop(
    "One or more final figure files are missing or too small:\n",
    paste(
      missing_figure_outputs,
      collapse = "\n"
    )
  )
}

cat(
  "\n========================================\n",
  "FIGURES 4-2 AND 4-3 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Figure 4-2 domain: Full Virginia\n",
  "Figure 4-3 domain: Northern Virginia zoom\n",
  "Figure 4-2 PNG method: ",
  statewide_png_method,
  "\nFigure 4-3 PNG method: ",
  nova_png_method,
  "\nResolution: ",
  png_dpi,
  " dpi\n\nFinal figure files:\n",
  paste(
    combined_figure_outputs,
    collapse = "\n"
  ),
  "\n\nShared QC/data directory:\n",
  out_dir,
  "\n",
  sep = ""
)
