#!/usr/bin/env Rscript

# ==============================================================================
# Figure 4-2 — final coordinate-anchored statistics boxes, lowered
# (a) Gasoline-attributable PM2.5, 2011-2020 mean
# (b) Diesel-attributable PM2.5, 2011-2020 mean
#
# CONUS mean, IQR (Q1-Q3), and maximum are placed in a wide white box inside
# the lower-right corner of each map. The box is anchored 0.40 degrees above
# the map's lower boundary, so it remains inside the panel while sitting below
# the mapped CONUS land area. Panel tags occupy a separate top strip.
# Source identities should be defined in the figure caption.
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(USAboundaries)
  library(scales)
  library(viridis)
})

# ------------------------------------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------------------------------------

input_dir <- "/scratch/xshan2/R_Code/Automobiles/FIGURES/NATIVE_CMAQ_REBUILT"

gasoline_file <- file.path(input_dir, "DECADE_ONR_mean_2011_2020.rds")
diesel_file   <- file.path(input_dir, "DECADE_NRD_mean_2011_2020.rds")

output_dir <- file.path(
  input_dir,
  "FIGURE_4_2_FINAL_PERFECT_INSIDE_BOX_LOWERED"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file <- file.path(
  output_dir,
  "Figure_4_2_FINAL_perfect_inside_box_lowered.pdf"
)
png_file <- file.path(
  output_dir,
  "Figure_4_2_FINAL_perfect_inside_box_lowered.png"
)
qc_file <- file.path(output_dir, "Figure_4_2_QC.csv")
plot_data_file <- file.path(output_dir, "Figure_4_2_plot_data.rds")

required_inputs <- c(gasoline_file, diesel_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

pm25_min <- 0
pm25_max <- 1.0
pm25_breaks <- seq(0, 1.0, by = 0.2)

conus_xlim <- c(-125.5, -64.5)
conus_ylim <- c(19.0, 50.5)

fig_width_in  <- 13.8
fig_height_in <- 5.80
png_dpi <- 600

# ------------------------------------------------------------------------------
# 2. Read the two decade-mean source surfaces
# ------------------------------------------------------------------------------

read_surface <- function(file_path, source_code, source_name) {
  dt <- as.data.table(readRDS(file_path))

  required_columns <- c("x", "y", "value")
  missing_columns <- setdiff(required_columns, names(dt))
  if (length(missing_columns) > 0) {
    stop(
      "Input file is missing required columns:\n",
      file_path,
      "\nMissing: ",
      paste(missing_columns, collapse = ", ")
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
    stop("No finite values were found in:\n", file_path)
  }
  if (anyDuplicated(out[, .(lon, lat)]) > 0) {
    stop("Duplicate longitude-latitude coordinates were found in:\n", file_path)
  }

  out[, source_code := source_code]
  out[, source_name := source_name]
  out
}

gasoline_dt <- read_surface(
  gasoline_file,
  source_code = "ONR",
  source_name = "Gasoline"
)
diesel_dt <- read_surface(
  diesel_file,
  source_code = "NRD",
  source_name = "Diesel"
)

if (!identical(
  gasoline_dt[order(lon, lat), .(lon, lat)],
  diesel_dt[order(lon, lat), .(lon, lat)]
)) {
  stop("Gasoline and diesel decade surfaces do not use identical coordinates.")
}

plot_data <- rbindlist(
  list(gasoline_dt, diesel_dt),
  use.names = TRUE
)

# ------------------------------------------------------------------------------
# 3. Retain grid-cell centers within CONUS
# ------------------------------------------------------------------------------

states <- USAboundaries::us_states(resolution = "low")
states <- states[
  !states$state_abbr %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS"),
]
states <- st_transform(states, 4326)

virginia <- states[states$state_abbr == "VA", ]
if (nrow(virginia) != 1) {
  stop("Virginia was not uniquely identified.")
}

grid_coordinates <- unique(plot_data[, .(lon, lat)])
setorder(grid_coordinates, lon, lat)

grid_points_ll <- st_as_sf(
  grid_coordinates,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)
grid_points_5070 <- st_transform(grid_points_ll, 5070)
states_5070 <- st_transform(states, 5070)

inside_conus <- lengths(st_intersects(grid_points_5070, states_5070)) > 0

conus_keys <- as.data.table(
  st_drop_geometry(grid_points_ll[inside_conus, ])
)[, .(lon, lat)]

setkey(plot_data, lon, lat)
setkey(conus_keys, lon, lat)
conus_data <- plot_data[conus_keys, nomatch = 0]

if (nrow(conus_data) == 0) {
  stop("No CMAQ grid-cell centers were retained for CONUS.")
}

# ------------------------------------------------------------------------------
# 4. Summary statistics and QC
# ------------------------------------------------------------------------------

qc <- conus_data[, .(
  n_grid_cells = .N,
  mean = mean(concentration, na.rm = TRUE),
  median = median(concentration, na.rm = TRUE),
  p05 = as.numeric(quantile(concentration, 0.05, na.rm = TRUE)),
  p25 = as.numeric(quantile(concentration, 0.25, na.rm = TRUE)),
  p75 = as.numeric(quantile(concentration, 0.75, na.rm = TRUE)),
  iqr = IQR(concentration, na.rm = TRUE, type = 7),
  p95 = as.numeric(quantile(concentration, 0.95, na.rm = TRUE)),
  p99 = as.numeric(quantile(concentration, 0.99, na.rm = TRUE)),
  maximum = max(concentration, na.rm = TRUE),
  percent_above_display_maximum =
    100 * mean(concentration > pm25_max, na.rm = TRUE)
), by = .(source_code, source_name)]

qc[, source_order := match(source_code, c("ONR", "NRD"))]
setorder(qc, source_order)
qc[, source_order := NULL]

fwrite(qc, qc_file)
saveRDS(
  list(
    conus = conus_data,
    summary_statistics = qc,
    display_limits = c(pm25_min, pm25_max),
    display_breaks = pm25_breaks
  ),
  plot_data_file
)

cat(
  "\n========================================\n",
  "FIGURE 4-2 QC SUMMARY\n",
  "========================================\n",
  sep = ""
)
print(qc)

# ------------------------------------------------------------------------------
# 5. Map panels
# ------------------------------------------------------------------------------

validate_summary_row <- function(summary_row) {
  if (nrow(summary_row) != 1) {
    stop("Each panel must receive exactly one summary-statistics row.")
  }
  invisible(TRUE)
}

make_map_panel <- function(data, summary_row = NULL, show_legend = FALSE) {
  p <- ggplot() +
    geom_raster(
      data = data,
      aes(x = lon, y = lat, fill = concentration),
      interpolate = FALSE
    ) +
    geom_sf(
      data = states,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey62",
      linewidth = 0.17
    ) +
    geom_sf(
      data = virginia,
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.55
    ) +
    coord_sf(
      xlim = conus_xlim,
      ylim = conus_ylim,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +
    scale_fill_viridis_c(
      option = "E",
      direction = 1,
      limits = c(pm25_min, pm25_max),
      breaks = pm25_breaks,
      labels = scales::label_number(accuracy = 0.1),
      oob = scales::squish,
      na.value = "white",
      name = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = grid::unit(82, "mm"),
        barheight = grid::unit(4.4, "mm"),
        frame.colour = "grey35",
        frame.linewidth = 0.30,
        ticks.colour = "grey35"
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 10, base_family = "sans") +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(
        fill = NA,
        color = "grey45",
        linewidth = 0.30
      ),
      legend.position = if (show_legend) "bottom" else "none",
      legend.justification = "center",
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size = 8.7),
      legend.margin = margin(t = 2, r = 0, b = 2, l = 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
    )

  if (!is.null(summary_row)) {
    validate_summary_row(summary_row)

    # All coordinates below are deliberately inside conus_xlim/conus_ylim.
    # This binds the statistics box to the actual sf panel instead of the
    # larger viewport that contains the aspect-ratio-constrained map.
    box_xmin <- -89.5
    box_xmax <- -66.2

    # Keep the box 0.40 degrees above the lower map border. Relative to the
    # preceding version, this moves the complete box down by 0.40 degrees
    # without changing its height or allowing it to cross the panel boundary.
    box_bottom_margin <- 0.40
    box_height <- 4.80
    box_ymin <- conus_ylim[1] + box_bottom_margin
    box_ymax <- box_ymin + box_height

    if (
      box_ymin <= conus_ylim[1] ||
        box_ymax >= conus_ylim[2]
    ) {
      stop("Statistics box coordinates fall outside the map panel.")
    }

    label_x <- -77.6
    value_x <- -67.2
    row_y <- box_ymax - c(0.85, 2.40, 3.95)

    label_text <- c("Mean", "IQR (Q1-Q3)", "Maximum")
    value_text <- c(
      sprintf("%.3f", summary_row$mean),
      sprintf("%.3f-%.3f", summary_row$p25, summary_row$p75),
      sprintf("%.3f", summary_row$maximum)
    )

    p <- p +
      annotate(
        "rect",
        xmin = box_xmin,
        xmax = box_xmax,
        ymin = box_ymin,
        ymax = box_ymax,
        fill = scales::alpha("white", 0.94),
        color = "grey30",
        linewidth = 0.28
      ) +
      annotate(
        "text",
        x = rep(label_x, 3),
        y = row_y,
        label = label_text,
        hjust = 1,
        vjust = 0.5,
        family = "sans",
        fontface = "plain",
        size = 3.28,
        color = "grey12"
      ) +
      annotate(
        "text",
        x = rep(value_x, 3),
        y = row_y,
        label = value_text,
        hjust = 1,
        vjust = 0.5,
        family = "sans",
        fontface = "plain",
        size = 3.28,
        color = "grey12"
      )
  }

  p
}

p_a <- make_map_panel(
  conus_data[source_code == "ONR"],
  summary_row = qc[source_code == "ONR"],
  show_legend = FALSE
)
p_b <- make_map_panel(
  conus_data[source_code == "NRD"],
  summary_row = qc[source_code == "NRD"],
  show_legend = FALSE
)
p_legend <- make_map_panel(
  conus_data[source_code == "ONR"],
  summary_row = NULL,
  show_legend = TRUE
) +
  theme(panel.border = element_blank())

get_legend <- function(plot_object) {
  g <- ggplotGrob(plot_object)
  guide_ids <- grep("^guide-box", g$layout$name)
  if (length(guide_ids) == 0) {
    stop("No legend guide-box was found.")
  }

  for (i in guide_ids) {
    candidate <- g$grobs[[i]]
    if (!inherits(candidate, "zeroGrob")) {
      return(candidate)
    }
  }

  stop("All legend guide-box objects were empty.")
}

legend_grob <- get_legend(p_legend)
g_a <- ggplotGrob(p_a)
g_b <- ggplotGrob(p_b)

# ------------------------------------------------------------------------------
# 6. Composite: panel tags, map panels, and shared legend
# ------------------------------------------------------------------------------

draw_tagged_panel <- function(grob_object, tag) {
  inner_layout <- grid::grid.layout(
    nrow = 2,
    ncol = 1,
    heights = grid::unit.c(
      grid::unit(0.31, "in"),
      grid::unit(1, "null")
    )
  )

  grid::pushViewport(grid::viewport(layout = inner_layout))

  grid::pushViewport(
    grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  grid::grid.text(
    label = tag,
    x = grid::unit(0.002, "npc"),
    y = grid::unit(0.44, "npc"),
    just = c("left", "center"),
    gp = grid::gpar(
      fontfamily = "sans",
      fontface = "bold",
      fontsize = 12.5
    )
  )
  grid::popViewport()

  grid::pushViewport(
    grid::viewport(layout.pos.row = 2, layout.pos.col = 1)
  )
  grid::grid.draw(grob_object)
  grid::popViewport(2)
}

draw_figure <- function() {
  grid::grid.newpage()

  page_layout <- grid::grid.layout(
    nrow = 3,
    ncol = 3,
    widths = grid::unit.c(
      grid::unit(0.16, "in"),
      grid::unit(1, "null"),
      grid::unit(0.16, "in")
    ),
    heights = grid::unit.c(
      grid::unit(0.16, "in"),
      grid::unit(1, "null"),
      grid::unit(0.12, "in")
    )
  )

  grid::pushViewport(grid::viewport(layout = page_layout))
  grid::pushViewport(
    grid::viewport(layout.pos.row = 2, layout.pos.col = 2)
  )

  outer_layout <- grid::grid.layout(
    nrow = 2,
    ncol = 3,
    widths = grid::unit.c(
      grid::unit(1, "null"),
      grid::unit(0.22, "in"),
      grid::unit(1, "null")
    ),
    heights = grid::unit.c(
      grid::unit(1, "null"),
      grid::unit(0.72, "in")
    )
  )

  grid::pushViewport(grid::viewport(layout = outer_layout))

  grid::pushViewport(
    grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
  )
  draw_tagged_panel(
    g_a,
    "(a)"
  )
  grid::popViewport()

  grid::pushViewport(
    grid::viewport(layout.pos.row = 1, layout.pos.col = 3)
  )
  draw_tagged_panel(
    g_b,
    "(b)"
  )
  grid::popViewport()

  grid::pushViewport(
    grid::viewport(layout.pos.row = 2, layout.pos.col = 1:3)
  )
  grid::grid.draw(legend_grob)
  grid::popViewport()

  grid::popViewport(3)
}

# ------------------------------------------------------------------------------
# 7. Save vector PDF
# ------------------------------------------------------------------------------

unlink(c(pdf_file, png_file), force = TRUE)

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
invisible(grDevices::dev.off())
Sys.sleep(1)

if (
  !file.exists(pdf_file) ||
    !is.finite(file.info(pdf_file)$size) ||
    file.info(pdf_file)$size <= 1000
) {
  stop("PDF creation failed:\n", pdf_file)
}

cat(
  "\nVector PDF created:\n",
  pdf_file,
  "\nSize: ",
  round(file.info(pdf_file)$size / 1024^2, 3),
  " MB\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 8. Create 600-dpi PNG
# ------------------------------------------------------------------------------

png_is_valid <- function(path) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size > 1000
}

png_method <- NA_character_

if (requireNamespace("ragg", quietly = TRUE)) {
  try({
    ragg::agg_png(
      filename = png_file,
      width = fig_width_in,
      height = fig_height_in,
      units = "in",
      res = png_dpi,
      background = "white"
    )
    draw_figure()
    invisible(grDevices::dev.off())
    Sys.sleep(1)
    if (png_is_valid(png_file)) {
      png_method <- "ragg"
    }
  }, silent = TRUE)
}

if (is.na(png_method) && isTRUE(capabilities("cairo"))) {
  if (file.exists(png_file)) {
    unlink(png_file, force = TRUE)
  }
  try({
    grDevices::png(
      filename = png_file,
      width = as.integer(round(fig_width_in * png_dpi)),
      height = as.integer(round(fig_height_in * png_dpi)),
      units = "px",
      res = png_dpi,
      type = "cairo-png",
      bg = "white",
      pointsize = 12
    )
    draw_figure()
    invisible(grDevices::dev.off())
    Sys.sleep(1)
    if (png_is_valid(png_file)) {
      png_method <- "base cairo-png"
    }
  }, silent = TRUE)
}

run_external <- function(command, args, method_name) {
  if (file.exists(png_file)) {
    unlink(png_file, force = TRUE)
  }

  output <- tryCatch(
    system2(
      command = command,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      structure(conditionMessage(e), status = 999L)
    }
  )

  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  Sys.sleep(1)

  if (status == 0L && png_is_valid(png_file)) {
    message("PNG method succeeded: ", method_name)
    return(TRUE)
  }

  message(
    method_name,
    " failed:\n",
    paste(output, collapse = "\n")
  )
  FALSE
}

find_existing_command <- function(command_name, fixed_paths = character()) {
  candidates <- unique(c(
    unname(Sys.which(command_name)),
    fixed_paths
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates[file.exists(candidates)][1]
}

if (is.na(png_method)) {
  pdftocairo_bin <- find_existing_command(
    "pdftocairo",
    c("/usr/bin/pdftocairo", "/bin/pdftocairo")
  )

  if (length(pdftocairo_bin) == 1 && !is.na(pdftocairo_bin)) {
    output_prefix <- tools::file_path_sans_ext(png_file)
    ok <- run_external(
      command = pdftocairo_bin,
      args = c(
        "-png",
        "-singlefile",
        "-r",
        as.character(png_dpi),
        pdf_file,
        output_prefix
      ),
      method_name = "pdftocairo"
    )
    if (ok) {
      png_method <- "pdftocairo"
    }
  }
}

if (is.na(png_method)) {
  pdftoppm_bin <- find_existing_command(
    "pdftoppm",
    c("/usr/bin/pdftoppm", "/bin/pdftoppm")
  )

  if (length(pdftoppm_bin) == 1 && !is.na(pdftoppm_bin)) {
    output_prefix <- tools::file_path_sans_ext(png_file)
    ok <- run_external(
      command = pdftoppm_bin,
      args = c(
        "-png",
        "-singlefile",
        "-r",
        as.character(png_dpi),
        pdf_file,
        output_prefix
      ),
      method_name = "pdftoppm"
    )
    if (ok) {
      png_method <- "pdftoppm"
    }
  }
}

if (is.na(png_method)) {
  stop(
    "The vector PDF was created successfully, but a PNG could not be created ",
    "in this Hopper session.\n\nPDF:\n",
    pdf_file
  )
}

# ------------------------------------------------------------------------------
# 9. Final verification
# ------------------------------------------------------------------------------

expected_outputs <- c(
  pdf_file,
  png_file,
  qc_file,
  plot_data_file
)

missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
empty_outputs <- expected_outputs[
  file.exists(expected_outputs) &
    (
      !is.finite(file.info(expected_outputs)$size) |
        file.info(expected_outputs)$size <= 0
    )
]

if (length(missing_outputs) > 0) {
  stop(
    "Missing expected outputs:\n",
    paste(missing_outputs, collapse = "\n")
  )
}
if (length(empty_outputs) > 0) {
  stop(
    "Empty expected outputs:\n",
    paste(empty_outputs, collapse = "\n")
  )
}

output_summary <- data.table(
  file = expected_outputs,
  size_mb = round(file.info(expected_outputs)$size / 1024^2, 3)
)

cat(
  "\n========================================\n",
  "FIGURE 4-2 COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "PNG method used: ",
  png_method,
  "\n",
  sep = ""
)
print(output_summary)
