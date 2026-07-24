#!/usr/bin/env Rscript

# ==============================================================================
# Figure 4-5: four-panel plot-only rebuild
# Fixed-2019 traffic proxy with P50 support regularization, Virginia 2011-2020
# ==============================================================================
#
# This script does not repeat the 2011-2020 downscaling calculation. It reads
# the validated decade-mean FST produced by:
#
#   Figure_4_5_Fixed_2019_AADT_P50_2011_2020.R
#
# Figure 4-4 already presents the wind-only fields. To avoid duplicating those
# maps, this version of Figure 4-5 contains only:
#
#   (a) Gasoline: wind + traffic
#   (b) Gasoline: traffic-induced change from wind only
#   (c) Diesel:   wind + traffic
#   (d) Diesel:   traffic-induced change from wind only
#
# Existing six-panel files are preserved. New output names end in "_4panel".
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(fst)
  library(sf)
  library(ggplot2)
  library(scales)
  library(grid)
  library(USAboundaries)
})

run_figure_4_5_four_panel <- function() {

  # ---------------------------------------------------------------------------
  # 1. Paths and fixed display settings
  # ---------------------------------------------------------------------------

  base_dir <- "/scratch/xshan2/R_Code/Automobiles"

  out_dir <- file.path(
    base_dir,
    "FIGURES",
    "NATIVE_CMAQ_REBUILT",
    "FIGURE_4_5_FIXED_2019_AADT_P50_FINAL_2011_2020"
  )

  decade_file <- file.path(
    out_dir,
    "Figure_4_5_P50_fixed2019_AADT_2011_2020_mean.fst"
  )

  previous_scale_file <- file.path(
    out_dir,
    "Figure_4_5_display_scales.csv"
  )

  pdf_file <- file.path(
    out_dir,
    "Figure_4_5_fixed2019_AADT_P50_VA_2011_2020_4panel.pdf"
  )

  png_file <- file.path(
    out_dir,
    "Figure_4_5_fixed2019_AADT_P50_VA_2011_2020_4panel.png"
  )

  tiff_file <- file.path(
    out_dir,
    "Figure_4_5_fixed2019_AADT_P50_VA_2011_2020_4panel.tiff"
  )

  scale_file <- file.path(
    out_dir,
    "Figure_4_5_4panel_display_scales.csv"
  )

  caption_file <- file.path(
    out_dir,
    "Figure_4_5_4panel_caption.txt"
  )

  completion_file <- file.path(
    out_dir,
    "Figure_4_5_4PANEL_PLOT_COMPLETE.txt"
  )

  figure_width_in <- 8.8
  figure_height_in <- 6.5
  raster_dpi <- 600

  absolute_upper_quantile <- 0.995
  difference_upper_quantile <- 0.995
  numerical_tolerance <- 1e-8

  p4s <- paste(
    "+proj=lcc", "+lat_1=33", "+lat_2=45", "+lat_0=40",
    "+lon_0=-97", "+a=6370000", "+b=6370000"
  )

  if (!file.exists(decade_file)) {
    stop(
      "Missing validated decade-mean input:\n",
      decade_file,
      "\nRun Figure_4_5_Fixed_2019_AADT_P50_2011_2020.R first."
    )
  }

  if (file.exists(completion_file)) {
    unlink(completion_file)
  }

  # ---------------------------------------------------------------------------
  # 2. Read and validate the completed decade-mean result
  # ---------------------------------------------------------------------------

  plot_data <- as.data.table(read_fst(decade_file))

  required_columns <- c(
    "grid_id",
    "x",
    "y",
    "parent_id",
    "source_code",
    "pm_wind",
    "pm_final",
    "delta_from_wind",
    "n_years"
  )

  missing_columns <- setdiff(required_columns, names(plot_data))

  if (length(missing_columns) > 0L) {
    stop(
      "The decade-mean input is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (
    nrow(plot_data) == 0L ||
      any(!is.finite(plot_data$x)) ||
      any(!is.finite(plot_data$y)) ||
      any(!is.finite(plot_data$pm_wind)) ||
      any(!is.finite(plot_data$pm_final)) ||
      any(!is.finite(plot_data$delta_from_wind))
  ) {
    stop("The decade-mean input is empty or contains non-finite plot values.")
  }

  if (!setequal(unique(plot_data$source_code), c("ONR", "NRD"))) {
    stop(
      "Expected source_code values ONR and NRD; found: ",
      paste(sort(unique(plot_data$source_code)), collapse = ", ")
    )
  }

  source_counts <- plot_data[, .(
    n_cells = .N,
    n_grid_id = uniqueN(grid_id),
    minimum_years = min(n_years),
    maximum_years = max(n_years)
  ), by = source_code]

  if (
    uniqueN(source_counts$n_cells) != 1L ||
      uniqueN(source_counts$n_grid_id) != 1L ||
      any(source_counts$n_cells != source_counts$n_grid_id) ||
      any(source_counts$minimum_years != 10L) ||
      any(source_counts$maximum_years != 10L)
  ) {
    stop("The two source fields do not contain the same complete 10-year grid.")
  }

  identity_error <- max(abs(
    plot_data$pm_final -
      plot_data$pm_wind -
      plot_data$delta_from_wind
  ))

  if (identity_error > numerical_tolerance) {
    stop(
      "The decade-mean delta identity failed. Maximum error = ",
      format(identity_error, scientific = TRUE)
    )
  }

  gasoline_data <- plot_data[source_code == "ONR"]
  diesel_data <- plot_data[source_code == "NRD"]

  # ---------------------------------------------------------------------------
  # 3. Virginia boundary and display limits
  # ---------------------------------------------------------------------------

  states_ll <- USAboundaries::us_states(resolution = "low")
  virginia_lcc <- states_ll[states_ll$state_abbr == "VA", ]
  virginia_lcc <- st_transform(st_as_sf(virginia_lcc), st_crs(p4s))

  if (nrow(virginia_lcc) != 1L) {
    stop("Could not resolve a single Virginia state boundary.")
  }

  va_bbox <- st_bbox(virginia_lcc)
  va_xlim <- c(va_bbox[["xmin"]], va_bbox[["xmax"]])
  va_ylim <- c(va_bbox[["ymin"]], va_bbox[["ymax"]])

  absolute_limit <- NA_real_
  difference_limit <- NA_real_
  scale_source <- "recalculated from validated decade-mean data"

  if (file.exists(previous_scale_file)) {
    previous_scales <- tryCatch(
      fread(previous_scale_file),
      error = function(e) NULL
    )

    if (
      !is.null(previous_scales) &&
        all(c("scale", "upper") %in% names(previous_scales)) &&
        any(previous_scales$scale == "absolute") &&
        any(previous_scales$scale == "difference")
    ) {
      absolute_limit <- previous_scales[
        scale == "absolute",
        as.numeric(upper[1])
      ]
      difference_limit <- abs(previous_scales[
        scale == "difference",
        as.numeric(upper[1])
      ])
      scale_source <- "reused from the validated six-panel figure"
    }
  }

  if (!is.finite(absolute_limit) || absolute_limit <= 0) {
    absolute_values <- c(plot_data$pm_wind, plot_data$pm_final)
    absolute_limit <- as.numeric(quantile(
      absolute_values,
      absolute_upper_quantile,
      names = FALSE,
      na.rm = TRUE
    ))
  }

  if (!is.finite(absolute_limit) || absolute_limit <= 0) {
    absolute_limit <- max(plot_data$pm_final, na.rm = TRUE)
  }

  if (!is.finite(difference_limit) || difference_limit <= 0) {
    difference_limit <- as.numeric(quantile(
      abs(plot_data$delta_from_wind),
      difference_upper_quantile,
      names = FALSE,
      na.rm = TRUE
    ))
  }

  if (!is.finite(difference_limit) || difference_limit <= 0) {
    difference_limit <- max(
      abs(plot_data$delta_from_wind),
      na.rm = TRUE
    )
  }

  if (
    !is.finite(absolute_limit) ||
      absolute_limit <= 0 ||
      !is.finite(difference_limit) ||
      difference_limit <= 0
  ) {
    stop("Could not derive valid concentration display limits.")
  }

  absolute_breaks <- pretty(c(0, absolute_limit), n = 5)
  absolute_breaks <- absolute_breaks[
    absolute_breaks >= 0 & absolute_breaks <= absolute_limit
  ]

  # Five explicitly symmetric breaks prevent rounded duplicate labels.
  difference_breaks <- seq(
    -difference_limit,
    difference_limit,
    length.out = 5L
  )

  fwrite(
    data.table(
      scale = c("absolute", "difference"),
      lower = c(0, -difference_limit),
      upper = c(absolute_limit, difference_limit),
      breaks = c(
        paste(signif(absolute_breaks, 8), collapse = "; "),
        paste(signif(difference_breaks, 8), collapse = "; ")
      ),
      label_accuracy = c(0.1, 0.01),
      scale_source = scale_source
    ),
    scale_file
  )

  # ---------------------------------------------------------------------------
  # 4. Publication map functions
  # ---------------------------------------------------------------------------

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
        linewidth = 0.30
      ),
      coord_sf(
        crs = st_crs(p4s),
        xlim = va_xlim,
        ylim = va_ylim,
        expand = FALSE,
        datum = NA,
        clip = "on"
      ),
      theme_void(base_size = 9, base_family = "sans"),
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(
          fill = NA,
          color = "grey45",
          linewidth = 0.25
        ),
        plot.title = element_text(
          size = 11.2,
          face = "bold",
          hjust = 0,
          margin = margin(b = 2.0)
        ),
        plot.margin = margin(t = 2, r = 4, b = 1, l = 4, unit = "pt")
      )
    )
  }

  make_absolute_panel <- function(
      panel_data,
      panel_title,
      show_legend = FALSE) {

    ggplot() +
      geom_raster(
        data = panel_data,
        aes(x = x, y = y, fill = pm_final),
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
          barwidth = unit(60, "mm"),
          barheight = unit(4.2, "mm"),
          frame.colour = "grey35",
          frame.linewidth = 0.25,
          ticks.colour = "grey35"
        )
      ) +
      labs(title = panel_title, x = NULL, y = NULL) +
      theme(
        legend.position = if (show_legend) "bottom" else "none",
        legend.title = element_text(size = 9.1),
        legend.text = element_text(size = 8.3)
      )
  }

  make_difference_panel <- function(
      panel_data,
      panel_title,
      show_legend = FALSE) {

    ggplot() +
      geom_raster(
        data = panel_data,
        aes(x = x, y = y, fill = delta_from_wind),
        interpolate = FALSE
      ) +
      common_map_layers() +
      scale_fill_gradientn(
        colours = difference_colours,
        values = rescale(c(
          -difference_limit,
          -difference_limit * 0.66,
          -difference_limit * 0.33,
          0,
          difference_limit * 0.33,
          difference_limit * 0.66,
          difference_limit
        )),
        limits = c(-difference_limit, difference_limit),
        breaks = difference_breaks,
        labels = label_number(accuracy = 0.01),
        oob = squish,
        na.value = "white",
        name = expression(Delta*PM[2.5]~(mu*g~m^{-3})),
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          title.hjust = 0.5,
          label.position = "bottom",
          barwidth = unit(60, "mm"),
          barheight = unit(4.2, "mm"),
          frame.colour = "grey35",
          frame.linewidth = 0.25,
          ticks.colour = "grey35"
        )
      ) +
      labs(title = panel_title, x = NULL, y = NULL) +
      theme(
        legend.position = if (show_legend) "bottom" else "none",
        legend.title = element_text(size = 9.1),
        legend.text = element_text(size = 8.3)
      )
  }

  panels <- list(
    make_absolute_panel(
      gasoline_data,
      "(a)"
    ),
    make_difference_panel(
      gasoline_data,
      "(b)"
    ),
    make_absolute_panel(
      diesel_data,
      "(c)"
    ),
    make_difference_panel(
      diesel_data,
      "(d)"
    )
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
    make_absolute_panel(gasoline_data, "", TRUE) +
      theme(
        panel.border = element_blank(),
        plot.title = element_blank()
      )
  )

  difference_legend <- get_legend(
    make_difference_panel(gasoline_data, "", TRUE) +
      theme(
        panel.border = element_blank(),
        plot.title = element_blank()
      )
  )

  panel_grobs <- lapply(panels, ggplotGrob)

  draw_main_figure <- function() {
    grid.newpage()

    figure_layout <- grid.layout(
      nrow = 3,
      ncol = 2,
      widths = unit(c(1, 1), "null"),
      heights = unit.c(
        unit(1, "null"),
        unit(1, "null"),
        unit(0.64, "in")
      )
    )

    pushViewport(viewport(layout = figure_layout))

    draw_at <- function(grob, row, column) {
      pushViewport(viewport(
        layout.pos.row = row,
        layout.pos.col = column
      ))
      grid.draw(grob)
      popViewport()
    }

    draw_at(panel_grobs[[1]], 1, 1)
    draw_at(panel_grobs[[2]], 1, 2)
    draw_at(panel_grobs[[3]], 2, 1)
    draw_at(panel_grobs[[4]], 2, 2)
    draw_at(absolute_legend, 3, 1)
    draw_at(difference_legend, 3, 2)

    popViewport()
  }

  # ---------------------------------------------------------------------------
  # 5. Save vector and 600-dpi raster outputs
  # ---------------------------------------------------------------------------

  pdf_partial <- paste0(pdf_file, ".partial")

  pdf(
    pdf_partial,
    width = figure_width_in,
    height = figure_height_in,
    useDingbats = FALSE,
    family = "sans"
  )
  draw_main_figure()
  dev.off()

  if (!file.exists(pdf_partial) || file.info(pdf_partial)$size < 10000) {
    stop("The vector figure was not rendered correctly: ", pdf_partial)
  }

  if (!file.copy(pdf_partial, pdf_file, overwrite = TRUE)) {
    stop("Could not install the completed vector figure: ", pdf_file)
  }
  unlink(pdf_partial)

  save_raster_figure <- function(file_path, file_type) {
    partial_path <- paste0(file_path, ".partial")
    ok <- FALSE

    if (file.exists(partial_path)) {
      unlink(partial_path)
    }

    if (requireNamespace("ragg", quietly = TRUE)) {
      tryCatch({
        if (file_type == "png") {
          ragg::agg_png(
            partial_path,
            width = figure_width_in,
            height = figure_height_in,
            units = "in",
            res = raster_dpi,
            background = "white"
          )
        } else {
          ragg::agg_tiff(
            partial_path,
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
        ok <- file.exists(partial_path) &&
          file.info(partial_path)$size > 0
      }, error = function(e) {
        warning(
          "ragg could not create ", file_type, ": ",
          conditionMessage(e)
        )
        while (dev.cur() > 1) dev.off()
      })
    }

    if (!ok && capabilities("cairo")) {
      if (file.exists(partial_path)) {
        unlink(partial_path)
      }

      tryCatch({
        if (file_type == "png") {
          png(
            partial_path,
            width = figure_width_in,
            height = figure_height_in,
            units = "in",
            res = raster_dpi,
            type = "cairo-png",
            bg = "white"
          )
        } else {
          tiff(
            partial_path,
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
        ok <- file.exists(partial_path) &&
          file.info(partial_path)$size > 0
      }, error = function(e) {
        warning(
          "Cairo could not create ", file_type, ": ",
          conditionMessage(e)
        )
        while (dev.cur() > 1) dev.off()
      })
    }

    if (ok) {
      if (!file.copy(partial_path, file_path, overwrite = TRUE)) {
        warning("Could not install completed ", file_type, ": ", file_path)
        ok <- FALSE
      }
    }

    if (file.exists(partial_path)) {
      unlink(partial_path)
    }

    ok
  }

  png_created <- save_raster_figure(png_file, "png")
  tiff_created <- save_raster_figure(tiff_file, "tiff")

  # ---------------------------------------------------------------------------
  # 6. Updated caption and completion record
  # ---------------------------------------------------------------------------

  caption_text <- paste0(
    "Figure 4-5 | Traffic-based refinement of wind-adjusted mobile-source ",
    "PM2.5 across Virginia, 2011-2020. Mean gasoline-attributable (a) and ",
    "diesel-attributable (c) PM2.5 after traffic adjustment, with ",
    "corresponding changes from the wind-only fields (b,d). Red and blue ",
    "indicate increases and decreases, respectively."
  )

  writeLines(caption_text, caption_file)

  completion_text <- c(
    "FIGURE 4-5 FOUR-PANEL PLOT COMPLETE",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input:", decade_file),
    paste("PDF:", pdf_file),
    paste(
      "PNG:",
      if (png_created) png_file else "not created; use the PDF"
    ),
    paste(
      "TIFF:",
      if (tiff_created) tiff_file else "not created; use the PDF"
    ),
    paste("Caption:", caption_file),
    paste("Scale record:", scale_file),
    paste(
      "Maximum decade identity error:",
      format(identity_error, scientific = TRUE)
    ),
    paste("Absolute display limit:", signif(absolute_limit, 8)),
    paste("Difference display limit:", signif(difference_limit, 8)),
    paste("Scale source:", scale_source),
    "The original six-panel outputs were not overwritten."
  )

  writeLines(completion_text, completion_file)

  cat(
    "\n============================================================\n",
    "FIGURE 4-5 FOUR-PANEL PLOT COMPLETE\n",
    "============================================================\n",
    paste(completion_text, collapse = "\n"),
    "\n",
    sep = ""
  )

  invisible(list(
    pdf_file = pdf_file,
    png_file = if (png_created) png_file else NA_character_,
    tiff_file = if (tiff_created) tiff_file else NA_character_,
    caption_file = caption_file,
    scale_file = scale_file,
    completion_file = completion_file
  ))
}

run_figure_4_5_four_panel()
