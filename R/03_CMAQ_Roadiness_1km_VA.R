#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fst)
  library(ggplot2)
  library(USAboundaries)
  library(viridis)
  library(scales)
})

# =========================================================
# SETTINGS
# =========================================================
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

base_dir   <- "/scratch/xshan2/R_Code/Automobiles"
annual_dir <- file.path(base_dir, "FIGURES")

roadiness_fst <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/VA/roadiness_1km_hw_loc_VA.fst"

years  <- 2011:2020
metric <- "mean"   # 如需 median：改成 "median"

# roadiness column mapping (固定保留 NRD/ONR/TOTAL 三层)
road_col_NRD <- "nroad.dist2_loc"
road_col_ONR <- "nroad.dist2_hw"

# outputs (only final useful data + 3 PDFs)
out_dir   <- file.path(base_dir, "DOWNSCALE_1KM_VA_FINAL")
out_data  <- file.path(out_dir, "DATA_FST")
out_fig   <- file.path(out_dir, "FIGURES")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)

# numeric safety (not truncation): avoid divide-by-zero
eps <- 1e-12

# Plot settings
facet_ncol <- 5
legend_upper_q <- 0.99  # only affects colorbar range (visualization), does NOT modify data

# =========================================================
# HELPERS
# =========================================================
read_annual_rds <- function(sector, metric, yr) {
  f <- file.path(annual_dir, sprintf("ANNUAL_%s_%s_%d.rds", sector, metric, yr))
  if (!file.exists(f)) stop("Missing annual RDS: ", f)
  dt <- readRDS(f)
  setDT(dt)
  if (!all(c("lon","lat","val") %in% names(dt))) stop("Annual RDS must have lon, lat, val: ", f)
  dt <- dt[is.finite(val) & !is.na(val)]
  dt
}

load_roadiness <- function(colname) {
  rd <- read_fst(roadiness_fst, as.data.table = TRUE)
  setDT(rd)
  if (!all(c("x","y") %in% names(rd))) stop("Roadiness fst must contain x,y.")
  if (!(colname %in% names(rd))) stop("Roadiness col not found: ", colname)

  out <- rd[, .(x, y, road = as.numeric(get(colname)))]
  out[!is.finite(road) | is.na(road), road := 0]
  out
}

make_pmgrid_and_mapidx <- function(nrd_dt, road_sf) {
  # Build unique CMAQ grid from NRD (assumed same grid for ONR/TOTAL)
  pmgrid <- unique(nrd_dt[, .(lon, lat)])
  pmgrid[, grid_id := .I]

  pm_sf <- st_as_sf(pmgrid, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pm_sf <- st_transform(pm_sf, st_crs(p4s))

  # For nearest_feature, both must be same CRS
  idx <- st_nearest_feature(road_sf, pm_sf)   # returns row index into pmgrid
  list(pmgrid = pmgrid, idx = idx)
}

downscale_one_layer <- function(pm_dt, pmgrid, idx, road_dt) {
  # pm_dt: lon, lat, val (pm12)
  # pmgrid: lon, lat, grid_id
  # idx: for each road cell, nearest pmgrid row
  # road_dt: x,y,road (already LCC), same order as used to create road_sf for idx

  # attach grid_id to pm values
  pm_dt2 <- merge(pm_dt[, .(lon, lat, pm12 = val)], pmgrid, by = c("lon","lat"), all.x = TRUE)
  if (any(is.na(pm_dt2$grid_id))) {
    stop("Some CMAQ cells could not be matched to pmgrid by lon/lat. Check grid consistency.")
  }
  setkey(pm_dt2, grid_id)

  # road -> grid_id via idx
  dt <- copy(road_dt)
  dt[, grid_id := pmgrid$grid_id[idx]]

  # attach pm12 by grid_id
  setkey(dt, grid_id)
  dt <- pm_dt2[dt]   # brings pm12

  # pure normalization weights (no truncation, no smoothing)
  dt[, road_eff := pmax(road, 0)]
  dt[, road_mean := mean(road_eff, na.rm = TRUE), by = grid_id]
  dt[, weight := (road_eff + eps) / (road_mean + eps)]
  dt[, pm1 := pm12 * weight]

  dt[, .(x, y, pm1)]
}

save_fst_xyval <- function(dt, outfile) {
  # dt: x,y,valcol
  write_fst(dt, outfile)
  cat("WROTE: ", outfile, " | rows=", nrow(dt), " | size=", file.info(outfile)$size, "\n", sep="")
}

read_stack_years <- function(layer = c("NRD","ONR","TOTAL"), years, metric) {
  layer <- match.arg(layer)

  files <- switch(layer,
    NRD   = file.path(out_data, sprintf("VA_1km_NRD_%s_%d.fst", metric, years)),
    ONR   = file.path(out_data, sprintf("VA_1km_ONR_%s_%d.fst", metric, years)),
    TOTAL = file.path(out_data, sprintf("VA_1km_TOTAL_%s_%d.fst", metric, years))
  )

  lst <- lapply(seq_along(years), function(i){
    f <- files[i]
    if (!file.exists(f)) stop("Missing: ", f)
    dt <- read_fst(f, as.data.table = TRUE)
    setDT(dt)
    if (layer == "TOTAL") setnames(dt, "pm1_total", "val")
    else setnames(dt, "pm1", "val")
    dt[, year := years[i]]
    dt
  })
  rbindlist(lst, use.names = TRUE)
}

plot_panel_va <- function(stack_dt, out_pdf, title_text = NULL) {
  # VA outline in p4s
  va_sf <- USAboundaries::us_states(resolution = "low")
  va_sf <- va_sf[va_sf$state_abbr == "VA", ]
  va_sf <- st_transform(st_as_sf(va_sf), st_crs(p4s))

  stack_dt <- stack_dt[is.finite(val) & !is.na(val)]
  stack_dt[, year_f := factor(year, levels = sort(unique(year)))]

  xlim <- range(stack_dt$x, na.rm = TRUE)
  ylim <- range(stack_dt$y, na.rm = TRUE)

  lim_upper <- as.numeric(quantile(stack_dt$val, legend_upper_q, na.rm = TRUE))
  if (!is.finite(lim_upper) || lim_upper <= 0) lim_upper <- max(stack_dt$val, na.rm = TRUE)

  p <- ggplot() +
    geom_raster(data = stack_dt, aes(x = x, y = y, fill = val)) +
    geom_sf(data = va_sf, fill = NA, color = "grey35", linewidth = 0.2) +
    coord_sf(crs = st_crs(p4s), xlim = xlim, ylim = ylim, expand = FALSE) +
    facet_wrap(~ year_f, ncol = facet_ncol) +
    scale_fill_viridis_c(
      name = expression(PM[2.5]~(mu*g/m^3)),
      limits = c(0, lim_upper),
      oob = scales::squish,
      na.value = "transparent"
    ) +
    # labs(title = title_text) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )

  grDevices::pdf(out_pdf, width = 16, height = 6.5, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  cat("WROTE: ", out_pdf, " | size=", file.info(out_pdf)$size, "\n", sep="")
}

# =========================================================
# RUN: produce ONLY final useful data (NRD/ONR/TOTAL) + 3 panel PDFs
# =========================================================

cat("=== VA 1km downscale (KEEP ONLY NRD/ONR/TOTAL FST; NO intermediate files) ===\n")
cat("Output root: ", out_dir, "\n", sep="")

# load roadiness once (two versions)
road_nrd <- load_roadiness(road_col_NRD)
road_onr <- load_roadiness(road_col_ONR)

# build sf points once for NN mapping (use NRD road grid coordinates as reference)
road_sf_ref <- st_as_sf(road_nrd[, .(x,y)], coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

for (yr in years) {
  cat("\n-----------------------------\n")
  cat("YEAR: ", yr, " | metric: ", metric, "\n", sep="")
  cat("-----------------------------\n")

  # Read annual CMAQ
  nrd12 <- read_annual_rds("NRD", metric, yr)
  onr12 <- read_annual_rds("ONR", metric, yr)

  # Build pm grid + NN mapping index (based on NRD grid lon/lat)
  mm <- make_pmgrid_and_mapidx(nrd12, road_sf_ref)
  pmgrid <- mm$pmgrid
  idx    <- mm$idx

  # Downscale NRD
  nrd1 <- downscale_one_layer(nrd12, pmgrid, idx, road_nrd)
  f_nrd <- file.path(out_data, sprintf("VA_1km_NRD_%s_%d.fst", metric, yr))
  save_fst_xyval(nrd1, f_nrd)

  # Downscale ONR (use same pmgrid + same idx)
  onr1 <- downscale_one_layer(onr12, pmgrid, idx, road_onr)
  f_onr <- file.path(out_data, sprintf("VA_1km_ONR_%s_%d.fst", metric, yr))
  save_fst_xyval(onr1, f_onr)

  # TOTAL = NRD + ONR (same x,y template)
  setkey(nrd1, x, y)
  setkey(onr1, x, y)
  tot1 <- nrd1[onr1]
  tot1[, pm1_total := pm1 + i.pm1]
  tot1 <- tot1[, .(x, y, pm1_total)]

  f_tot <- file.path(out_data, sprintf("VA_1km_TOTAL_%s_%d.fst", metric, yr))
  save_fst_xyval(tot1, f_tot)

  # minimal console diagnostics (no files)
  cat("Ranges (NRD pm1):   ", paste(range(nrd1$pm1, na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("Ranges (ONR pm1):   ", paste(range(onr1$pm1, na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("Ranges (TOTAL pm1): ", paste(range(tot1$pm1_total, na.rm=TRUE), collapse=" ~ "), "\n", sep="")

  rm(nrd12, onr12, mm, pmgrid, idx, nrd1, onr1, tot1); gc()
}

cat("\n=== MAKING PANEL MAPS (2011–2020) ===\n")

stack_nrd <- read_stack_years("NRD", years, metric)
stack_onr <- read_stack_years("ONR", years, metric)
stack_tot <- read_stack_years("TOTAL", years, metric)

pdf_nrd <- file.path(out_fig, sprintf("PANEL_VA_NRD_%s_2011_2020.pdf", metric))
pdf_onr <- file.path(out_fig, sprintf("PANEL_VA_ONR_%s_2011_2020.pdf", metric))
pdf_tot <- file.path(out_fig, sprintf("PANEL_VA_TOTAL_%s_2011_2020.pdf", metric))

plot_panel_va(stack_nrd, pdf_nrd, title_text = NULL)
plot_panel_va(stack_onr, pdf_onr, title_text = NULL)
plot_panel_va(stack_tot, pdf_tot, title_text = NULL)

cat("\nALL DONE ✅\n")
cat("Final data (FST): ", out_data, "\n", sep="")
cat("Panel maps (PDF): ", out_fig,  "\n", sep="")
