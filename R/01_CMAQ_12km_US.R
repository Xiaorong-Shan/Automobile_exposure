# ====================================================================================================
# Original data: /home/xshan2/HAQ_LAB/Sumaiya/cmaq/cmaq_output/R/data/conus_data/with_filter/
# Auto PM2.5 (CMAQ) — Annual Mean/Median (2011–2020)
# Inputs (daily CSVs):
#   /scratch/xshan2/R_Code/Automobiles/NRD/PM25_CONUS_NRD_YYYY.csv
#   /scratch/xshan2/R_Code/Automobiles/ONR/PM25_CONUS_ONR_YYYY.csv
# Outputs:
#   1) Annual grids (RDS):
#      FIGURES/ANNUAL_NRD_mean_YYYY.rds
#      FIGURES/ANNUAL_NRD_median_YYYY.rds
#      FIGURES/ANNUAL_ONR_mean_YYYY.rds
#      FIGURES/ANNUAL_ONR_median_YYYY.rds
#      FIGURES/ANNUAL_TOTAL_mean_YYYY.rds
#      FIGURES/ANNUAL_TOTAL_median_YYYY.rds
#   2) Bad-day logs (CSV):
#      FIGURES/BAD_DAYS_YYYY.csv  (union bad days for NRD/ONR)
#   3) Facet panel maps (PDF):
#      FIGURES/PANEL_MAPS/PANEL_<SECTOR>_<metric>_2011_2020.pdf
#   4) Boxplots (PDF):
#      FIGURES/BOXPLOTS/BOXPLOT_<metric>_2011_2020.pdf
# Notes:
#   - "Bad days" removed using conservative thresholds on daily p99/max.
#   - Uses geom_raster (fast) + USAboundaries state outlines (WGS84).
#   - All plot titles are commented out (per request).
# ====================================================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(USAboundaries)
  library(scales)
  library(viridis)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
base_dir <- "/scratch/xshan2/R_Code/Automobiles"
in_dir_nrd <- file.path(base_dir, "NRD")
in_dir_onr <- file.path(base_dir, "ONR")

years <- 2011:2020

out_annual_dir <- file.path(base_dir, "FIGURES")
out_panel_dir  <- file.path(base_dir, "FIGURES", "PANEL_MAPS")
out_box_dir    <- file.path(base_dir, "FIGURES", "BOXPLOTS")
dir.create(out_annual_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_panel_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_box_dir,    showWarnings = FALSE, recursive = TRUE)

# WGS84 extent for CONUS plotting
xlim_conus <- c(-125, -66)
ylim_conus <- c(24, 50)

# Bad-day detection thresholds (tuned for CMAQ blow-ups)
max_cut <- 200
p99_cut <- 50

# Sampling for daily QC (bigger => more accurate; smaller => faster)
sample_n <- 2e6

# Limits in maps: use global upper quantile across 2011–2020 for each plot
legend_upper_q <- 0.99

# For boxplots: sample cells per year/sector/metric to avoid huge memory
box_sample_per_year <- 200000

# -----------------------------
# BASE MAP: CONUS STATES (WGS84)
# -----------------------------
states_sf <- USAboundaries::us_states(resolution = "low")
states_sf <- states_sf[!states_sf$state_abbr %in% c("AK","HI","PR","VI","GU","MP","AS"), ]
states_sf <- st_transform(states_sf, 4326)

# -----------------------------
# HELPERS
# -----------------------------
read_sector_year_csv <- function(sector = c("NRD","ONR"), yr) {
  sector <- match.arg(sector)
  f <- if (sector == "NRD") {
    file.path(in_dir_nrd, sprintf("PM25_CONUS_NRD_%d.csv", yr))
  } else {
    file.path(in_dir_onr, sprintf("PM25_CONUS_ONR_%d.csv", yr))
  }
  if (!file.exists(f)) stop("Missing file: ", f)

  dt <- fread(
    f,
    na.strings = c("NA", "", "NaN", "Inf", "-Inf"),
    showProgress = FALSE
  )

  # Expected columns: lon, lat, Date, value (and maybe sector)
  # Normalize column names
  if (!("value" %in% names(dt))) stop("Expected column 'value' not found in: ", f)
  if (!("lon" %in% names(dt)))   stop("Expected column 'lon' not found in: ", f)
  if (!("lat" %in% names(dt)))   stop("Expected column 'lat' not found in: ", f)
  if (!("Date" %in% names(dt)))  stop("Expected column 'Date' not found in: ", f)

  dt[, Date := as.IDate(Date)]
  dt[, value := as.numeric(value)]
  dt <- dt[is.finite(value) & !is.na(value)]

  # Physical cleanup
  dt[value < 0, value := 0]

  dt[, sector := sector]
  dt[, year := yr]
  dt
}

daily_stats_sample <- function(dt, sample_n = 2e6) {
  # dt needs: Date, value
  setDT(dt)
  dt[, {
    v <- value
    v <- v[is.finite(v)]
    if (length(v) > sample_n) v <- sample(v, sample_n)
    if (length(v) == 0) {
      list(p99 = NA_real_, mx = NA_real_)
    } else {
      list(
        p99 = as.numeric(quantile(v, 0.99, na.rm = TRUE)),
        mx  = as.numeric(max(v, na.rm = TRUE))
      )
    }
  }, by = Date]
}

compute_annual_mean_median <- function(dt_day, bad_dates_union = NULL) {
  # dt_day columns: lon, lat, Date, value
  setDT(dt_day)
  if (!is.null(bad_dates_union) && length(bad_dates_union) > 0) {
    dt_day <- dt_day[!Date %in% bad_dates_union]
  }
  # Annual by grid cell
  dt_mean <- dt_day[, .(val = mean(value, na.rm = TRUE)), by = .(lon, lat)]
  dt_med  <- dt_day[, .(val = median(value, na.rm = TRUE)), by = .(lon, lat)]
  list(mean = dt_mean, median = dt_med)
}

save_annual_rds <- function(dt_grid, sector, metric, yr) {
  f <- file.path(out_annual_dir, sprintf("ANNUAL_%s_%s_%d.rds", sector, metric, yr))
  saveRDS(dt_grid[, .(lon, lat, val)], f)
  f
}

write_bad_days_log <- function(yr, bad_nrd, bad_onr, bad_union) {
  f <- file.path(out_annual_dir, sprintf("BAD_DAYS_%d.csv", yr))
  out <- data.table(
    year = yr,
    Date = sort(unique(as.IDate(bad_union))),
    bad_in = "NRD_or_ONR"
  )
  fwrite(out, f)
  # Also print short summary
  cat(sprintf("\nBad day summary %d:\n", yr))
  cat("  NRD bad days: ", if (length(bad_nrd)) paste(bad_nrd, collapse=", ") else "(none)", "\n", sep="")
  cat("  ONR bad days: ", if (length(bad_onr)) paste(bad_onr, collapse=", ") else "(none)", "\n", sep="")
  cat("  UNION bad days:", if (length(bad_union)) paste(bad_union, collapse=", ") else "(none)", "\n", sep="")
  cat("  Wrote log: ", f, "\n", sep="")
  f
}

# Stack annual RDS for panel maps / boxplots
read_annual_stack <- function(sector = c("NRD","ONR","TOTAL"),
                              metric = c("mean","median"),
                              years = 2011:2020) {
  sector <- match.arg(sector)
  metric <- match.arg(metric)
  lst <- lapply(years, function(yr){
    f <- file.path(out_annual_dir, sprintf("ANNUAL_%s_%s_%d.rds", sector, metric, yr))
    if (!file.exists(f)) stop("Missing annual RDS: ", f)
    dt <- readRDS(f)
    setDT(dt)
    dt[, year := yr]
    dt[, sector := sector]
    dt[, metric := metric]
    dt[, .(lon, lat, val, year, sector, metric)]
  })
  rbindlist(lst, use.names = TRUE)
}

plot_panel_map <- function(stack_dt, out_pdf, ncol = 5) {
  lim_upper <- as.numeric(quantile(stack_dt$val, legend_upper_q, na.rm = TRUE))
  if (!is.finite(lim_upper) || lim_upper <= 0) lim_upper <- max(stack_dt$val, na.rm = TRUE)

  stack_dt[, year_f := factor(year, levels = sort(unique(year)))]

  p <- ggplot() +
    geom_raster(data = stack_dt, aes(x = lon, y = lat, fill = val)) +
    geom_sf(data = states_sf, fill = NA, color = "grey35", linewidth = 0.15) +
    coord_sf(xlim = xlim_conus, ylim = ylim_conus, expand = FALSE) +
    facet_wrap(~ year_f, ncol = ncol) +
    scale_fill_viridis_c(
      name   = expression(PM[2.5]~(mu*g/m^3)),
      limits = c(0, lim_upper),
      oob    = scales::squish,
      na.value = "transparent"
    ) +
    # labs(title = "COMMENTED OUT") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      # plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )

  grDevices::pdf(out_pdf, width = 16, height = 6.5, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  cat("WROTE: ", out_pdf, " | size=", file.info(out_pdf)$size, "\n", sep="")
}

plot_boxplot <- function(dt_vals, out_pdf) {
  # dt_vals columns: year, sector, metric, val
  dt_vals[, year_f := factor(year, levels = sort(unique(year)))]
  dt_vals[, sector := factor(sector, levels = c("NRD","ONR","TOTAL"))]

  # Use a log-like view only if needed; keep linear here
  p <- ggplot(dt_vals, aes(x = year_f, y = val)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.10, width = 0.6) +
    facet_wrap(~ sector, ncol = 1, scales = "free_y") +
    # labs(title = "COMMENTED OUT") +
    labs(x = NULL, y = expression(PM[2.5]~(mu*g/m^3))) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      # plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  grDevices::pdf(out_pdf, width = 12, height = 10, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  cat("WROTE: ", out_pdf, " | size=", file.info(out_pdf)$size, "\n", sep="")
}

# -----------------------------
# 1) BUILD ANNUAL GRIDS (NRD/ONR/TOTAL) WITH BAD-DAY REMOVAL
# -----------------------------
for (yr in years) {
  cat("\n========================================\n")
  cat("YEAR: ", yr, "\n", sep="")
  cat("========================================\n")

  # Read daily CSVs
  nrd <- read_sector_year_csv("NRD", yr)
  onr <- read_sector_year_csv("ONR", yr)

  # Identify bad days (sector-specific)
  stat_nrd <- daily_stats_sample(nrd, sample_n = sample_n)
  stat_onr <- daily_stats_sample(onr, sample_n = sample_n)

  bad_nrd <- stat_nrd[(mx > max_cut) | (p99 > p99_cut), Date]
  bad_onr <- stat_onr[(mx > max_cut) | (p99 > p99_cut), Date]
  bad_union <- sort(unique(c(bad_nrd, bad_onr)))

  # Write bad-day log
  write_bad_days_log(yr, bad_nrd, bad_onr, bad_union)

  # Annual for each sector after removing union bad days
  ann_nrd <- compute_annual_mean_median(nrd, bad_union)
  ann_onr <- compute_annual_mean_median(onr, bad_union)

  # Save NRD/ONR annual grids
  f1 <- save_annual_rds(ann_nrd$mean,   "NRD", "mean",   yr)
  f2 <- save_annual_rds(ann_nrd$median, "NRD", "median", yr)
  f3 <- save_annual_rds(ann_onr$mean,   "ONR", "mean",   yr)
  f4 <- save_annual_rds(ann_onr$median, "ONR", "median", yr)
  cat("Saved:\n  ", f1, "\n  ", f2, "\n  ", f3, "\n  ", f4, "\n", sep="")

  # TOTAL = NRD + ONR on the same grid (after same bad days removed)
  setkey(ann_nrd$mean, lon, lat); setkey(ann_onr$mean, lon, lat)
  setkey(ann_nrd$median, lon, lat); setkey(ann_onr$median, lon, lat)

  tot_mean <- ann_nrd$mean[ann_onr$mean]
  tot_mean[, val := pmax(0, val) + pmax(0, i.val)]
  tot_mean <- tot_mean[, .(lon, lat, val)]

  tot_med <- ann_nrd$median[ann_onr$median]
  tot_med[, val := pmax(0, val) + pmax(0, i.val)]
  tot_med <- tot_med[, .(lon, lat, val)]

  f5 <- save_annual_rds(tot_mean, "TOTAL", "mean", yr)
  f6 <- save_annual_rds(tot_med,  "TOTAL", "median", yr)
  cat("Saved:\n  ", f5, "\n  ", f6, "\n", sep="")

  rm(nrd, onr, stat_nrd, stat_onr, ann_nrd, ann_onr, tot_mean, tot_med)
  gc()
}

cat("\nDONE: Annual grids written to:\n  ", out_annual_dir, "\n", sep="")

# -----------------------------
# 2) PANEL MAPS (2011–2020 in ONE figure per sector/metric)
# -----------------------------
make_panel <- function(sector, metric) {
  dt <- read_annual_stack(sector, metric, years)
  out_pdf <- file.path(out_panel_dir, sprintf("PANEL_%s_%s_2011_2020.pdf", sector, metric))
  plot_panel_map(dt, out_pdf, ncol = 5)
}

make_panel("NRD",   "mean")
make_panel("NRD",   "median")
make_panel("ONR",   "mean")
make_panel("ONR",   "median")
make_panel("TOTAL", "mean")
make_panel("TOTAL", "median")

# -----------------------------
# 3) BOXPLOTS (one PDF for mean, one for median; 2011–2020)
# -----------------------------
build_box_dt <- function(metric = c("mean","median")) {
  metric <- match.arg(metric)
  lst <- lapply(c("NRD","ONR","TOTAL"), function(sec){
    dt <- read_annual_stack(sec, metric, years)
    # Sample to keep Hopper safe
    dt <- dt[, .SD[sample.int(.N, size = min(.N, box_sample_per_year))], by = year]
    dt
  })
  rbindlist(lst, use.names = TRUE)
}

box_mean <- build_box_dt("mean")
box_med  <- build_box_dt("median")

plot_boxplot(box_mean, file.path(out_box_dir, "BOXPLOT_mean_2011_2020.pdf"))
plot_boxplot(box_med,  file.path(out_box_dir, "BOXPLOT_median_2011_2020.pdf"))

cat("\nALL DONE ✅\n",
    "Annual grids: ", out_annual_dir, "\n",
    "Panel maps:   ", out_panel_dir,  "\n",
    "Boxplots:     ", out_box_dir,    "\n", sep="")
