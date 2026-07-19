# ==============================================================================
# Figure 4-2 (FINAL, Hopper-safe)
# Native-resolution CMAQ-ISAM PM2.5 from gasoline and diesel vehicles
# 2011-2020 mean, shown for CONUS and Virginia
#
# Confirmed source definitions:
#   ONR = gasoline vehicles
#   NRD = diesel vehicles
#
# Figure panels (the figure itself shows panel tags only):
#   (a) Gasoline vehicles, CONUS
#   (b) Diesel vehicles, CONUS
#   (c) Gasoline vehicles, Virginia
#   (d) Diesel vehicles, Virginia
#
# Source/region descriptions belong in the figure caption, not inside the maps.
#
# This script intentionally does NOT require:
#   - package Cairo
#   - package lwgeom
#   - X11
#   - pdftoppm at a fixed path
#
# It always writes a vector PDF with base R. It then attempts several independent
# headless PNG routes in sequence and accepts the first one that succeeds:
#   1. base-R cairo PNG, when R was compiled with cairo support
#   2. ragg::agg_png, when ragg is installed
#   3. Cairo::CairoPNG, when Cairo is installed
#   4. grDevices::bitmap (Ghostscript-backed)
#   5. Ghostscript conversion from PDF
#   6. pdftocairo conversion from PDF
#   7. pdftoppm conversion from PDF
#   8. ImageMagick (magick/convert) conversion from PDF
#
# The script reports success ONLY after all required outputs actually exist.
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

input_dir <- paste0(
  "/scratch/xshan2/R_Code/Automobiles/FIGURES/",
  "NATIVE_CMAQ_REBUILT"
)

onr_file <- file.path(
  input_dir,
  "DECADE_ONR_mean_2011_2020.rds"
)

nrd_file <- file.path(
  input_dir,
  "DECADE_NRD_mean_2011_2020.rds"
)

out_dir <- file.path(
  input_dir,
  "FIGURE_4_2_FINAL"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf_file <- file.path(
  out_dir,
  "Figure_4_2_native_CMAQ_gasoline_diesel_CONUS_VA_ABCD_2011_2020.pdf"
)

png_file <- file.path(
  out_dir,
  "Figure_4_2_native_CMAQ_gasoline_diesel_CONUS_VA_ABCD_2011_2020.png"
)

qc_file <- file.path(
  out_dir,
  "Figure_4_2_QC.csv"
)

data_file <- file.path(
  out_dir,
  "Figure_4_2_plot_data.rds"
)

required_inputs <- c(onr_file, nrd_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing decade-mean input files:\n",
    paste(missing_inputs, collapse = "\n")
  )
}

# Figure dimensions. These proportions give two map rows plus one compact legend.
fig_width_in  <- 11.6
fig_height_in <- 6.25
png_dpi       <- 600

# Display values above this pooled CONUS quantile at the upper color limit.
# Numerical values are never altered.
upper_quantile <- 0.995

# Map windows
conus_xlim <- c(-125.0, -66.4)
conus_ylim <- c(24.0, 50.0)

# ------------------------------------------------------------------------------
# 2. Read source-specific decade means
# ------------------------------------------------------------------------------

read_decade_surface <- function(file_path, sector_code, source_name) {

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

  out[, sector_code := sector_code]
  out[, source_name := source_name]
  out
}

plot_data <- rbindlist(
  list(
    read_decade_surface(
      file_path = onr_file,
      sector_code = "ONR",
      source_name = "Gasoline vehicles"
    ),
    read_decade_surface(
      file_path = nrd_file,
      sector_code = "NRD",
      source_name = "Diesel vehicles"
    )
  ),
  use.names = TRUE
)

plot_data[, source_name := factor(
  source_name,
  levels = c("Gasoline vehicles", "Diesel vehicles")
)]

# Confirm that the two source surfaces use identical coordinates.
coord_onr <- plot_data[sector_code == "ONR", .(lon, lat)][order(lon, lat)]
coord_nrd <- plot_data[sector_code == "NRD", .(lon, lat)][order(lon, lat)]

if (!identical(coord_onr, coord_nrd)) {
  stop("ONR and NRD decade surfaces do not use identical grid coordinates.")
}

# ------------------------------------------------------------------------------
# 3. State boundaries and geographic filtering
# ------------------------------------------------------------------------------

states <- USAboundaries::us_states(resolution = "low")

states <- states[
  !states$state_abbr %in%
    c("AK", "HI", "PR", "VI", "GU", "MP", "AS"),
]

states <- st_transform(states, 4326)
va <- states[states$state_abbr == "VA", ]

if (nrow(va) != 1) {
  stop("Virginia was not uniquely identified in the state boundary data.")
}

# Perform point-in-polygon filtering in an equal-area projected CRS. This avoids
# longitude-latitude planar-operation warnings and does not require lwgeom.
grid_coords <- unique(plot_data[, .(lon, lat)])
setorder(grid_coords, lon, lat)

grid_points_ll <- st_as_sf(
  grid_coords,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

grid_points_5070 <- st_transform(grid_points_ll, 5070)
states_5070 <- st_transform(states, 5070)
va_5070 <- st_transform(va, 5070)

inside_conus <- lengths(
  st_intersects(grid_points_5070, states_5070)
) > 0

inside_va <- lengths(
  st_intersects(grid_points_5070, va_5070)
) > 0

conus_keys <- as.data.table(
  st_drop_geometry(grid_points_ll[inside_conus, ])
)[, .(lon, lat)]

va_keys <- as.data.table(
  st_drop_geometry(grid_points_ll[inside_va, ])
)[, .(lon, lat)]

setkey(plot_data, lon, lat)
setkey(conus_keys, lon, lat)
setkey(va_keys, lon, lat)

conus <- plot_data[conus_keys, nomatch = 0]
va_data <- plot_data[va_keys, nomatch = 0]

if (nrow(conus) == 0) {
  stop("No CMAQ grid-cell centers were retained for CONUS.")
}

if (nrow(va_data) == 0) {
  stop("No CMAQ grid-cell centers were retained for Virginia.")
}

# Virginia display window
va_bbox <- st_bbox(va)
va_xlim <- c(
  as.numeric(va_bbox["xmin"]) - 0.18,
  as.numeric(va_bbox["xmax"]) + 0.18
)
va_ylim <- c(
  as.numeric(va_bbox["ymin"]) - 0.13,
  as.numeric(va_bbox["ymax"]) + 0.13
)

# ------------------------------------------------------------------------------
# 4. Shared display scale and QC
# ------------------------------------------------------------------------------

# Use the pooled CONUS gasoline + diesel distribution so all four panels are
# directly comparable without allowing Virginia duplication to affect the scale.
upper_limit <- as.numeric(
  quantile(
    conus$concentration,
    probs = upper_quantile,
    na.rm = TRUE
  )
)

if (!is.finite(upper_limit) || upper_limit <= 0) {
  upper_limit <- max(conus$concentration, na.rm = TRUE)
}

legend_breaks <- pretty(c(0, upper_limit), n = 5)
legend_breaks <- legend_breaks[
  legend_breaks >= 0 & legend_breaks <= upper_limit
]

if (length(legend_breaks) < 3) {
  legend_breaks <- seq(0, upper_limit, length.out = 5)
}

region_stats <- function(dt, region_name) {

  ans <- dt[, .(
    n_grid_cells = .N,
    mean = mean(concentration, na.rm = TRUE),
    median = median(concentration, na.rm = TRUE),
    p05 = as.numeric(quantile(concentration, 0.05, na.rm = TRUE)),
    p25 = as.numeric(quantile(concentration, 0.25, na.rm = TRUE)),
    p75 = as.numeric(quantile(concentration, 0.75, na.rm = TRUE)),
    p95 = as.numeric(quantile(concentration, 0.95, na.rm = TRUE)),
    p99 = as.numeric(quantile(concentration, 0.99, na.rm = TRUE)),
    maximum = max(concentration, na.rm = TRUE),
    percent_above_display_limit = 100 * mean(
      concentration > upper_limit,
      na.rm = TRUE
    ),
    minimum_years = if (all(is.na(n_years))) {
      NA_integer_
    } else {
      min(n_years, na.rm = TRUE)
    },
    maximum_years = if (all(is.na(n_years))) {
      NA_integer_
    } else {
      max(n_years, na.rm = TRUE)
    }
  ), by = .(sector_code, source_name)]

  ans[, region := region_name]
  ans[, display_upper_quantile := upper_quantile]
  ans[, display_upper_limit := upper_limit]

  setcolorder(
    ans,
    c(
      "region",
      "sector_code",
      "source_name",
      setdiff(
        names(ans),
        c("region", "sector_code", "source_name")
      )
    )
  )

  ans
}

qc <- rbindlist(
  list(
    region_stats(conus, "CONUS"),
    region_stats(va_data, "Virginia")
  ),
  use.names = TRUE
)

fwrite(qc, qc_file)

saveRDS(
  list(
    conus = conus,
    virginia = va_data,
    display_upper_quantile = upper_quantile,
    display_upper_limit = upper_limit
  ),
  data_file
)

cat("\n========================================\n")
cat("FIGURE 4-2 DATA QC\n")
cat("========================================\n")
print(qc)
cat(
  "\nShared display upper limit: ",
  signif(upper_limit, 5),
  " ug/m3\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 5. Publication-style panels
# ------------------------------------------------------------------------------

make_panel <- function(
    data,
    boundary,
    xlim,
    ylim,
    highlight_va = FALSE,
    show_legend = FALSE) {

  p <- ggplot() +
    geom_raster(
      data = data,
      aes(
        x = lon,
        y = lat,
        fill = concentration
      ),
      interpolate = FALSE
    ) +
    geom_sf(
      data = boundary,
      inherit.aes = FALSE,
      fill = NA,
      color = "grey62",
      linewidth = 0.17
    )

  if (highlight_va) {
    p <- p + geom_sf(
      data = va,
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 0.60
    )
  }

  p +
    coord_sf(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE,
      datum = NA,
      clip = "on"
    ) +
    scale_fill_viridis_c(
      option = "E",
      direction = 1,
      limits = c(0, upper_limit),
      breaks = legend_breaks,
      labels = scales::label_number(accuracy = 0.1),
      oob = scales::squish,
      na.value = "white",
      name = expression(PM[2.5] ~ (mu * g ~ m^{-3})),
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom",
        barwidth = grid::unit(88, "mm"),
        barheight = grid::unit(4.2, "mm"),
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
      legend.position = if (show_legend) "bottom" else "none",
      legend.justification = "center",
      legend.title = element_text(size = 9.3),
      legend.text = element_text(size = 8.5),
      legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(1, 5, 3, 5, unit = "pt")
    )
}

p_a <- make_panel(
  data = conus[sector_code == "ONR"],
  boundary = states,
  xlim = conus_xlim,
  ylim = conus_ylim,
  highlight_va = TRUE,
  show_legend = FALSE
)

p_b <- make_panel(
  data = conus[sector_code == "NRD"],
  boundary = states,
  xlim = conus_xlim,
  ylim = conus_ylim,
  highlight_va = TRUE,
  show_legend = FALSE
)

# Neighboring-state outlines remain visible for geographic context, while the
# concentration fields in panels (c) and (d) are restricted to Virginia.
p_c <- make_panel(
  data = va_data[sector_code == "ONR"],
  boundary = states,
  xlim = va_xlim,
  ylim = va_ylim,
  highlight_va = TRUE,
  show_legend = FALSE
)

p_d <- make_panel(
  data = va_data[sector_code == "NRD"],
  boundary = states,
  xlim = va_xlim,
  ylim = va_ylim,
  highlight_va = TRUE,
  show_legend = FALSE
)

# Create one shared legend.
p_legend <- make_panel(
  data = conus[sector_code == "ONR"],
  boundary = states,
  xlim = conus_xlim,
  ylim = conus_ylim,
  highlight_va = FALSE,
  show_legend = TRUE
) +
  theme(
    panel.border = element_blank()
  )

get_legend <- function(plot_object) {

  grob <- ggplotGrob(plot_object)
  guide_ids <- grep("^guide-box", grob$layout$name)

  if (length(guide_ids) == 0) {
    stop("No legend guide-box was found.")
  }

  for (i in guide_ids) {
    candidate <- grob$grobs[[i]]
    if (!inherits(candidate, "zeroGrob")) {
      return(candidate)
    }
  }

  stop("All legend guide-box objects were empty.")
}

legend_grob <- get_legend(p_legend)

g_a <- ggplotGrob(p_a)
g_b <- ggplotGrob(p_b)
g_c <- ggplotGrob(p_c)
g_d <- ggplotGrob(p_d)

# ------------------------------------------------------------------------------
# 6. Composite figure drawing function
# ------------------------------------------------------------------------------

draw_figure <- function() {

  grid::grid.newpage()

  layout <- grid::grid.layout(
    nrow = 3,
    ncol = 2,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit.c(
      grid::unit(1.00, "null"),
      grid::unit(0.86, "null"),
      grid::unit(0.64, "in")
    )
  )

  grid::pushViewport(
    grid::viewport(layout = layout)
  )

  # Each panel cell contains a small tag row and the map itself. This keeps the
  # panel identifiers outside the mapped data and avoids all internal subtitles.
  draw_panel_with_tag <- function(grob_object, tag, row, col) {

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = row,
        layout.pos.col = col
      )
    )

    inner_layout <- grid::grid.layout(
      nrow = 2,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(0.19, "in"),
        grid::unit(1, "null")
      )
    )

    grid::pushViewport(
      grid::viewport(layout = inner_layout)
    )

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    grid::grid.text(
      label = tag,
      x = grid::unit(1.5, "mm"),
      y = grid::unit(0.5, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(
        fontfamily = "sans",
        fontface = "bold",
        fontsize = 11
      )
    )

    grid::popViewport()

    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
      )
    )

    grid::grid.draw(grob_object)

    grid::popViewport(3)
  }

  draw_panel_with_tag(g_a, "(a)", 1, 1)
  draw_panel_with_tag(g_b, "(b)", 1, 2)
  draw_panel_with_tag(g_c, "(c)", 2, 1)
  draw_panel_with_tag(g_d, "(d)", 2, 2)

  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 3,
      layout.pos.col = 1:2
    )
  )

  grid::grid.draw(legend_grob)

  grid::popViewport()
  grid::popViewport()
}

# ------------------------------------------------------------------------------
# 7. Save a vector PDF with base R (headless and already proven on Hopper)
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

if (!file.exists(pdf_file) || file.info(pdf_file)$size <= 1000) {
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
# 8. Hopper-safe PNG creation with multiple independent fallbacks
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

  if (file.exists(png_file)) {
    unlink(png_file, force = TRUE)
  }

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

  if (isTRUE(result)) {
    message("PNG method succeeded: ", method_name)
  }

  isTRUE(result)
}

png_method <- NA_character_

# Method 1: base-R cairo device, when compiled into this R installation.
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

# Method 2: ragg, when installed.
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

# Method 3: Cairo package, when installed. The script does not require it.
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

# Method 4: base-R bitmap. This may find Ghostscript through R_GSCMD even when
# Sys.which('gs') is empty.
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

# External conversion helper.
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
      structure(
        conditionMessage(e),
        status = 999L
      )
    }
  )

  status <- attr(output, "status")
  if (is.null(status)) status <- 0L

  if (status != 0L) {
    message(
      method_name,
      " failed with status ",
      status,
      ":\n",
      paste(output, collapse = "\n")
    )
  }

  Sys.sleep(1)

  if (status == 0L && png_is_valid(png_file)) {
    message("PNG method succeeded: ", method_name)
    return(TRUE)
  }

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

# Method 5: Ghostscript PDF -> PNG.
if (is.na(png_method)) {

  gs_bin <- find_existing_command(
    "gs",
    c("/usr/bin/gs", "/bin/gs", Sys.getenv("R_GSCMD"))
  )

  if (length(gs_bin) == 1 && !is.na(gs_bin)) {
    ok <- run_external(
      command = gs_bin,
      args = c(
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=png16m",
        paste0("-r", png_dpi),
        "-dTextAlphaBits=4",
        "-dGraphicsAlphaBits=4",
        paste0("-sOutputFile=", png_file),
        pdf_file
      ),
      method_name = "Ghostscript"
    )
    if (ok) png_method <- "Ghostscript"
  }
}

# Method 6: pdftocairo.
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
    if (ok) png_method <- "pdftocairo"
  }
}

# Method 7: pdftoppm.
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
        "-r",
        as.character(png_dpi),
        "-singlefile",
        pdf_file,
        output_prefix
      ),
      method_name = "pdftoppm"
    )
    if (ok) png_method <- "pdftoppm"
  }
}

# Method 8: ImageMagick already on PATH.
if (is.na(png_method)) {

  magick_bin <- find_existing_command(
    "magick",
    c("/usr/bin/magick", "/bin/magick")
  )

  if (length(magick_bin) == 1 && !is.na(magick_bin)) {
    ok <- run_external(
      command = magick_bin,
      args = c(
        "-density",
        as.character(png_dpi),
        paste0(pdf_file, "[0]"),
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
    if (ok) png_method <- "ImageMagick magick"
  }
}

if (is.na(png_method)) {

  convert_bin <- find_existing_command(
    "convert",
    c("/usr/bin/convert", "/bin/convert")
  )

  if (length(convert_bin) == 1 && !is.na(convert_bin)) {
    ok <- run_external(
      command = convert_bin,
      args = c(
        "-density",
        as.character(png_dpi),
        paste0(pdf_file, "[0]"),
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
    if (ok) png_method <- "ImageMagick convert"
  }
}

# Last fallback: attempt ImageMagick through a login shell/module environment.
if (is.na(png_method) && file.exists("/bin/bash")) {

  shell_command <- paste0(
    "source /etc/profile >/dev/null 2>&1 || true; ",
    "source /etc/profile.d/modules.sh >/dev/null 2>&1 || true; ",
    "module load imagemagick/7.1.0 >/dev/null 2>&1 || ",
    "module load imagemagick >/dev/null 2>&1 || ",
    "module load ImageMagick >/dev/null 2>&1 || true; ",
    "if command -v magick >/dev/null 2>&1; then ",
    "magick -density ", png_dpi, " '", pdf_file, "[0]' ",
    "-background white -alpha remove -alpha off '", png_file, "'; ",
    "elif command -v convert >/dev/null 2>&1; then ",
    "convert -density ", png_dpi, " '", pdf_file, "[0]' ",
    "-background white -alpha remove -alpha off '", png_file, "'; ",
    "else exit 127; fi"
  )

  output <- tryCatch(
    system2(
      "/bin/bash",
      args = c("-lc", shQuote(shell_command)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      structure(conditionMessage(e), status = 999L)
    }
  )

  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  Sys.sleep(1)

  if (status == 0L && png_is_valid(png_file)) {
    png_method <- "ImageMagick via module/login shell"
    message("PNG method succeeded: ", png_method)
  } else {
    message(
      "ImageMagick module fallback failed:\n",
      paste(output, collapse = "\n")
    )
  }
}

if (is.na(png_method) || !png_is_valid(png_file)) {
  stop(
    "The vector PDF was created successfully, but no available headless PNG ",
    "renderer/converter was found in this Hopper session.\n\n",
    "PDF:\n", pdf_file, "\n\n",
    "The script tried base cairo, ragg, Cairo, bitmap/Ghostscript, gs, ",
    "pdftocairo, pdftoppm, and ImageMagick.\n",
    "Run the following in the Hopper terminal to identify an available converter:\n",
    "command -v gs; command -v pdftocairo; command -v pdftoppm; ",
    "command -v magick; command -v convert\n"
  )
}

# ------------------------------------------------------------------------------
# 9. Final verification — no false success message
# ------------------------------------------------------------------------------

expected <- c(
  pdf_file,
  png_file,
  qc_file,
  data_file
)

missing <- expected[!file.exists(expected)]
empty <- expected[
  file.exists(expected) &
    (!is.finite(file.info(expected)$size) | file.info(expected)$size <= 0)
]

if (length(missing) > 0) {
  stop(
    "Missing expected outputs:\n",
    paste(missing, collapse = "\n")
  )
}

if (length(empty) > 0) {
  stop(
    "Empty expected outputs:\n",
    paste(empty, collapse = "\n")
  )
}

summary_out <- data.table(
  file = expected,
  size_mb = round(file.info(expected)$size / 1024^2, 3)
)

cat("\n========================================\n")
cat("FIGURE 4-2 COMPLETED SUCCESSFULLY\n")
cat("========================================\n")
cat("PNG method used: ", png_method, "\n", sep = "")
print(summary_out)
cat("\nOutput directory:\n", out_dir, "\n", sep = "")
