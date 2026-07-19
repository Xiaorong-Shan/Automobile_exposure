#!/usr/bin/env Rscript

# ==============================================================================
# Figure 4-2 — final compact two-panel version
# Native-resolution CMAQ-ISAM PM2.5 attributable to gasoline and diesel vehicles
# averaged over 2011–2020.
#
# Confirmed source definitions:
#   ONR = gasoline vehicles
#   NRD = diesel vehicles
#
# Panels:
#   (a) Gasoline-attributable PM2.5 across CONUS
#   (b) Diesel-attributable PM2.5 across CONUS
#
# Design revisions requested:
#   - Panel tags (a) and (b) are OUTSIDE the map boxes, directly above them.
#   - Shared legend has enough space so numbers are not clipped.
#   - Uses the earlier blue-grey-yellow palette (viridis option "E").
#   - Shared PM2.5 legend fixed at 0–1.0 ug m-3.
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

output_dir <- file.path(input_dir, "FIGURE_4_2_FINAL_AB_LAYOUTFIX")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pdf_file       <- file.path(output_dir, "Figure_4_2_native_CMAQ_gasoline_diesel_CONUS_AB_layoutfix_2011_2020.pdf")
png_file       <- file.path(output_dir, "Figure_4_2_native_CMAQ_gasoline_diesel_CONUS_AB_layoutfix_2011_2020.png")
qc_file        <- file.path(output_dir, "Figure_4_2_QC.csv")
plot_data_file <- file.path(output_dir, "Figure_4_2_plot_data.rds")

required_inputs <- c(gasoline_file, diesel_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input files:\n", paste(missing_inputs, collapse = "\n"))
}

pm25_min <- 0
pm25_max <- 1.0
pm25_breaks <- seq(0, 1.0, by = 0.2)

conus_xlim <- c(-125.0, -66.4)
conus_ylim <- c(24.0, 50.0)

fig_width_in  <- 11.5
fig_height_in <- 4.05
png_dpi <- 600

# ------------------------------------------------------------------------------
# 2. Read source-specific decade means
# ------------------------------------------------------------------------------

read_surface <- function(file_path, source_code, source_name) {
  dt <- as.data.table(readRDS(file_path))

  required_columns <- c("x", "y", "value")
  missing_columns <- setdiff(required_columns, names(dt))
  if (length(missing_columns) > 0) {
    stop("Input file is missing required columns:\n",
         file_path,
         "\nMissing: ",
         paste(missing_columns, collapse = ", "))
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

  out <- out[is.finite(lon) & is.finite(lat) & is.finite(concentration)]

  if (nrow(out) == 0) stop("No finite values were found in:\n", file_path)
  if (anyDuplicated(out[, .(lon, lat)]) > 0) {
    stop("Duplicate longitude-latitude coordinates were found in:\n", file_path)
  }

  out[, source_code := source_code]
  out[, source_name := source_name]
  out
}

gasoline_dt <- read_surface(gasoline_file, "ONR", "Gasoline vehicles")
diesel_dt   <- read_surface(diesel_file,   "NRD", "Diesel vehicles")

if (!identical(
  gasoline_dt[order(lon, lat), .(lon, lat)],
  diesel_dt[order(lon, lat),   .(lon, lat)]
)) {
  stop("Gasoline and diesel decade surfaces do not use identical coordinates.")
}

plot_data <- rbindlist(list(gasoline_dt, diesel_dt), use.names = TRUE)

# ------------------------------------------------------------------------------
# 3. CONUS and Virginia boundaries
# ------------------------------------------------------------------------------

states <- USAboundaries::us_states(resolution = "low")
states <- states[!states$state_abbr %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS"), ]
states <- st_transform(states, 4326)

virginia <- states[states$state_abbr == "VA", ]
if (nrow(virginia) != 1) stop("Virginia was not uniquely identified.")

# Filter grid centers to CONUS in a projected CRS
grid_coordinates <- unique(plot_data[, .(lon, lat)])
setorder(grid_coordinates, lon, lat)

grid_points_ll   <- st_as_sf(grid_coordinates, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
grid_points_5070 <- st_transform(grid_points_ll, 5070)
states_5070      <- st_transform(states, 5070)

inside_conus <- lengths(st_intersects(grid_points_5070, states_5070)) > 0

conus_keys <- as.data.table(st_drop_geometry(grid_points_ll[inside_conus, ]))[, .(lon, lat)]
setkey(plot_data, lon, lat)
setkey(conus_keys, lon, lat)
conus_data <- plot_data[conus_keys, nomatch = 0]

if (nrow(conus_data) == 0) stop("No CMAQ grid-cell centers were retained for CONUS.")

# ------------------------------------------------------------------------------
# 4. QC and saved plotting data
# ------------------------------------------------------------------------------

qc <- conus_data[, .(
  n_grid_cells = .N,
  mean    = mean(concentration, na.rm = TRUE),
  median  = median(concentration, na.rm = TRUE),
  p05     = as.numeric(quantile(concentration, 0.05, na.rm = TRUE)),
  p25     = as.numeric(quantile(concentration, 0.25, na.rm = TRUE)),
  p75     = as.numeric(quantile(concentration, 0.75, na.rm = TRUE)),
  p95     = as.numeric(quantile(concentration, 0.95, na.rm = TRUE)),
  p99     = as.numeric(quantile(concentration, 0.99, na.rm = TRUE)),
  maximum = max(concentration, na.rm = TRUE),
  percent_above_display_maximum = 100 * mean(concentration > pm25_max, na.rm = TRUE)
), by = .(source_code, source_name)]

fwrite(qc, qc_file)
saveRDS(list(conus = conus_data,
             display_limits = c(pm25_min, pm25_max),
             display_breaks = pm25_breaks),
        plot_data_file)

cat("\n========================================\n",
    "FIGURE 4-2 QC SUMMARY\n",
    "========================================\n", sep = "")
print(qc)

# ------------------------------------------------------------------------------
# 5. Map panels
# ------------------------------------------------------------------------------

make_panel <- function(data, show_legend = FALSE) {
  ggplot() +
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
      plot.background  = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(fill = NA, color = "grey45", linewidth = 0.30),
      legend.position  = if (show_legend) "bottom" else "none",
      legend.justification = "center",
      legend.title = element_text(size = 9.5),
      legend.text  = element_text(size = 8.7),
      legend.margin = margin(t = 2, r = 0, b = 2, l = 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 1, r = 4, b = 1, l = 4, unit = "pt")
    )
}

p_a <- make_panel(conus_data[source_code == "ONR"], show_legend = FALSE)
p_b <- make_panel(conus_data[source_code == "NRD"], show_legend = FALSE)
p_legend <- make_panel(conus_data[source_code == "ONR"], show_legend = TRUE) +
  theme(panel.border = element_blank())

get_legend <- function(plot_object) {
  g <- ggplotGrob(plot_object)
  guide_ids <- grep("^guide-box", g$layout$name)
  if (length(guide_ids) == 0) stop("No legend guide-box was found.")
  for (i in guide_ids) {
    candidate <- g$grobs[[i]]
    if (!inherits(candidate, "zeroGrob")) return(candidate)
  }
  stop("All legend guide-box objects were empty.")
}

legend_grob <- get_legend(p_legend)
g_a <- ggplotGrob(p_a)
g_b <- ggplotGrob(p_b)

# ------------------------------------------------------------------------------
# 6. Draw composite with panel tags OUTSIDE map boxes
# ------------------------------------------------------------------------------

draw_figure <- function() {
  grid::grid.newpage()

  outer_layout <- grid::grid.layout(
    nrow = 2,
    ncol = 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit.c(
      grid::unit(1, "null"),
      grid::unit(0.74, "in")
    )
  )

  grid::pushViewport(grid::viewport(layout = outer_layout))

  draw_panel_with_tag <- function(grob_object, tag, row, column) {
    grid::pushViewport(
      grid::viewport(layout.pos.row = row, layout.pos.col = column)
    )

    inner_layout <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(0.13, "in"),
        grid::unit(1, "null")
      )
    )

    grid::pushViewport(grid::viewport(layout = inner_layout))

    # Tag row OUTSIDE the panel box
    grid::pushViewport(
      grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
    )
    grid::grid.text(
      label = tag,
      x = grid::unit(0.01, "npc"),
      y = grid::unit(0.52, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = 11)
    )
    grid::popViewport()

    # Actual panel box
    grid::pushViewport(
      grid::viewport(layout.pos.row = 2, layout.pos.col = 1)
    )
    grid::grid.draw(grob_object)
    grid::popViewport(3)
  }

  draw_panel_with_tag(g_a, "(a)", 1, 1)
  draw_panel_with_tag(g_b, "(b)", 1, 2)

  grid::pushViewport(
    grid::viewport(layout.pos.row = 2, layout.pos.col = 1:2)
  )
  grid::grid.draw(legend_grob)
  grid::popViewport()

  grid::popViewport()
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

if (!file.exists(pdf_file) ||
    !is.finite(file.info(pdf_file)$size) ||
    file.info(pdf_file)$size <= 1000) {
  stop("PDF creation failed:\n", pdf_file)
}

cat("\nVector PDF created:\n",
    pdf_file,
    "\nSize: ",
    round(file.info(pdf_file)$size / 1024^2, 3),
    " MB\n", sep = "")

# ------------------------------------------------------------------------------
# 8. Hopper-safe PNG creation
# ------------------------------------------------------------------------------

png_is_valid <- function(path) {
  file.exists(path) &&
    is.finite(file.info(path)$size) &&
    file.info(path)$size > 1000
}

close_extra_devices <- function(start_device) {
  while (grDevices::dev.cur() > 1 && grDevices::dev.cur() != start_device) {
    try(grDevices::dev.off(), silent = TRUE)
  }
}

try_png_device <- function(method_name, open_device) {
  if (file.exists(png_file)) unlink(png_file, force = TRUE)
  start_device <- grDevices::dev.cur()

  result <- tryCatch({
    open_device()
    draw_figure()
    invisible(grDevices::dev.off())
    Sys.sleep(1)
    png_is_valid(png_file)
  }, error = function(e) {
    close_extra_devices(start_device)
    message(method_name, " failed: ", conditionMessage(e))
    FALSE
  })

  if (isTRUE(result)) message("PNG method succeeded: ", method_name)
  isTRUE(result)
}

png_method <- NA_character_

if (isTRUE(capabilities("cairo")) && is.na(png_method)) {
  ok <- try_png_device(
    "base grDevices::png(type='cairo-png')",
    function() {
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
    }
  )
  if (ok) png_method <- "base cairo-png"
}

if (requireNamespace("ragg", quietly = TRUE) && is.na(png_method)) {
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
  if (ok) png_method <- "ragg"
}

if (requireNamespace("Cairo", quietly = TRUE) && is.na(png_method)) {
  ok <- try_png_device(
    "Cairo::CairoPNG",
    function() {
      Cairo::CairoPNG(
        filename = png_file,
        width = as.integer(round(fig_width_in * png_dpi)),
        height = as.integer(round(fig_height_in * png_dpi)),
        unit = "px",
        pointsize = 12,
        bg = "white",
        res = png_dpi
      )
    }
  )
  if (ok) png_method <- "Cairo package"
}

if (is.na(png_method)) {
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
  if (ok) png_method <- "base bitmap / Ghostscript"
}

run_external <- function(command, args, method_name) {
  if (file.exists(png_file)) unlink(png_file, force = TRUE)

  output <- tryCatch(
    system2(command = command, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 999L)
  )

  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  if (status != 0L) {
    message(method_name, " failed with status ", status, ":\n", paste(output, collapse = "\n"))
  }

  Sys.sleep(1)

  if (status == 0L && png_is_valid(png_file)) {
    message("PNG method succeeded: ", method_name)
    return(TRUE)
  }
  FALSE
}

find_existing_command <- function(command_name, fixed_paths = character()) {
  candidates <- unique(c(unname(Sys.which(command_name)), fixed_paths))
  candidates <- candidates[nzchar(candidates)]
  candidates[file.exists(candidates)][1]
}

if (is.na(png_method)) {
  gs_bin <- find_existing_command("gs", c("/usr/bin/gs", "/bin/gs", Sys.getenv("R_GSCMD")))
  if (length(gs_bin) == 1 && !is.na(gs_bin)) {
    ok <- run_external(
      command = gs_bin,
      args = c("-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m",
               paste0("-r", png_dpi), "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
               paste0("-sOutputFile=", png_file), pdf_file),
      method_name = "Ghostscript"
    )
    if (ok) png_method <- "Ghostscript"
  }
}

if (is.na(png_method)) {
  pdftocairo_bin <- find_existing_command("pdftocairo", c("/usr/bin/pdftocairo", "/bin/pdftocairo"))
  if (length(pdftocairo_bin) == 1 && !is.na(pdftocairo_bin)) {
    output_prefix <- tools::file_path_sans_ext(png_file)
    ok <- run_external(
      command = pdftocairo_bin,
      args = c("-png", "-singlefile", "-r", as.character(png_dpi), pdf_file, output_prefix),
      method_name = "pdftocairo"
    )
    if (ok) png_method <- "pdftocairo"
  }
}

if (is.na(png_method)) {
  pdftoppm_bin <- find_existing_command("pdftoppm", c("/usr/bin/pdftoppm", "/bin/pdftoppm"))
  if (length(pdftoppm_bin) == 1 && !is.na(pdftoppm_bin)) {
    output_prefix <- tools::file_path_sans_ext(png_file)
    ok <- run_external(
      command = pdftoppm_bin,
      args = c("-png", "-r", as.character(png_dpi), "-singlefile", pdf_file, output_prefix),
      method_name = "pdftoppm"
    )
    if (ok) png_method <- "pdftoppm"
  }
}

if (is.na(png_method)) {
  magick_bin <- find_existing_command("magick", c("/usr/bin/magick", "/bin/magick"))
  if (length(magick_bin) == 1 && !is.na(magick_bin)) {
    ok <- run_external(
      command = magick_bin,
      args = c("-density", as.character(png_dpi), paste0(pdf_file, "[0]"),
               "-background", "white", "-alpha", "remove", "-alpha", "off", png_file),
      method_name = "ImageMagick magick"
    )
    if (ok) png_method <- "ImageMagick magick"
  }
}

if (is.na(png_method)) {
  convert_bin <- find_existing_command("convert", c("/usr/bin/convert", "/bin/convert"))
  if (length(convert_bin) == 1 && !is.na(convert_bin)) {
    ok <- run_external(
      command = convert_bin,
      args = c("-density", as.character(png_dpi), paste0(pdf_file, "[0]"),
               "-background", "white", "-alpha", "remove", "-alpha", "off", png_file),
      method_name = "ImageMagick convert"
    )
    if (ok) png_method <- "ImageMagick convert"
  }
}

if (is.na(png_method) && file.exists("/bin/bash")) {
  shell_command <- paste0(
    "source /etc/profile >/dev/null 2>&1 || true; ",
    "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; ",
    "module load imagemagick/7.1.0 >/dev/null 2>&1 || ",
    "module load imagemagick >/dev/null 2>&1 || ",
    "module load ImageMagick >/dev/null 2>&1 || true; ",
    "if command -v magick >/dev/null 2>&1; then ",
    "magick -density ", png_dpi, " '", pdf_file, "[0]' -background white -alpha remove -alpha off '", png_file, "'; ",
    "elif command -v convert >/dev/null 2>&1; then ",
    "convert -density ", png_dpi, " '", pdf_file, "[0]' -background white -alpha remove -alpha off '", png_file, "'; ",
    "else exit 127; fi"
  )

  shell_output <- tryCatch(
    system2("/bin/bash", args = c("-lc", shQuote(shell_command)), stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), status = 999L)
  )

  shell_status <- attr(shell_output, "status")
  if (is.null(shell_status)) shell_status <- 0L
  Sys.sleep(1)

  if (shell_status == 0L && png_is_valid(png_file)) {
    png_method <- "ImageMagick via module/login shell"
    message("PNG method succeeded: ", png_method)
  } else {
    message("ImageMagick module fallback failed:\n", paste(shell_output, collapse = "\n"))
  }
}

if (is.na(png_method) || !png_is_valid(png_file)) {
  stop("The vector PDF was created successfully, but no available headless PNG renderer/converter was found in this Hopper session.\n\nPDF:\n", pdf_file)
}

# ------------------------------------------------------------------------------
# 9. Final verification
# ------------------------------------------------------------------------------

expected_outputs <- c(pdf_file, png_file, qc_file, plot_data_file)

missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
empty_outputs <- expected_outputs[
  file.exists(expected_outputs) &
    (!is.finite(file.info(expected_outputs)$size) | file.info(expected_outputs)$size <= 0)
]

if (length(missing_outputs) > 0) {
  stop("Missing expected outputs:\n", paste(missing_outputs, collapse = "\n"))
}

if (length(empty_outputs) > 0) {
  stop("Empty expected outputs:\n", paste(empty_outputs, collapse = "\n"))
}

output_summary <- data.table(
  file = expected_outputs,
  size_mb = round(file.info(expected_outputs)$size / 1024^2, 3)
)

cat("\n========================================\n",
    "FIGURE 4-2 COMPLETED SUCCESSFULLY\n",
    "========================================\n",
    "PNG method used: ", png_method, "\n", sep = "")

print(output_summary)
cat("\nOutput directory:\n", output_dir, "\n", sep = "")
