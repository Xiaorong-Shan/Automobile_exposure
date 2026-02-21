#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fst)
  library(ggplot2)
  library(USAboundaries)
  library(viridis)
  library(scales)
  library(raster)
})

# =========================================================
# SETTINGS
# =========================================================
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

base_dir   <- "/scratch/xshan2/R_Code/Automobiles"
annual_dir <- file.path(base_dir, "FIGURES")

roadiness_fst <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/VA/roadiness_1km_hw_loc_VA.fst"

# Wind data folder and filename template
wind_dir <- "/scratch/xshan2/NLDAS_wind/VA"
wind_tif_fmt <- file.path(wind_dir, "WindDerived_VA_%d.tif")  # <-- adjust if needed

years  <- 2011:2020
metric <- "mean"

# roadiness column mapping
road_col_NRD <- "nroad.dist2_loc"
road_col_ONR <- "nroad.dist2_hw"

# meteo settings (ONLY wind speed term for now)
beta_wind  <- 1.0      # meteo = (u_ref / u)^beta
use_u_ref  <- TRUE     # TRUE: stabilize magnitude by u_ref (median); FALSE: pure 1/u^beta

# outputs
out_dir   <- file.path(base_dir, "DOWNSCALE_1KM_VA_METEO_RAW_NORM")
out_data  <- file.path(out_dir, "DATA_FST")
out_fig   <- file.path(out_dir, "FIGURES")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)

eps <- 1e-12

# Plot settings
facet_ncol <- 5
legend_upper_q <- 0.99

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

make_pmgrid_and_mapidx <- function(pm_dt, road_sf_ref) {
  pmgrid <- unique(pm_dt[, .(lon, lat)])
  pmgrid[, grid_id := .I]

  pm_sf <- st_as_sf(pmgrid, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pm_sf <- st_transform(pm_sf, st_crs(p4s))

  idx <- st_nearest_feature(road_sf_ref, pm_sf)   # length = nrow(road grid)
  list(pmgrid = pmgrid, idx = idx)
}

# wind brick: lon/lat WGS84; from your summary:
# layer2: WindDir_rad
# layer3: WindSpeed_eff (min=0.5)
extract_wind_to_road <- function(road_dt, yr, p4s, eps = 1e-12) {
  wind_tif <- sprintf(wind_tif_fmt, yr)
  if (!file.exists(wind_tif)) stop("Missing wind tif: ", wind_tif)

  wd <- raster::brick(wind_tif)

  # LCC (x,y) -> lon/lat for extraction
  pts_lcc <- st_as_sf(road_dt[, .(x,y)], coords = c("x","y"), crs = st_crs(p4s))
  pts_ll  <- st_transform(pts_lcc, 4326)
  ll      <- st_coordinates(pts_ll)
  lonlat  <- cbind(ll[,1], ll[,2])  # (lon, lat)

  wdir <- raster::extract(wd[[2]], lonlat)
  wspd <- raster::extract(wd[[3]], lonlat)

  out <- copy(road_dt)
  out[, `:=`(
    wdir = as.numeric(wdir),
    wspd = as.numeric(wspd),
    inv_u = 1 / (as.numeric(wspd) + eps)
  )]
  out
}

downscale_one_layer_road_meteo <- function(pm_dt, pmgrid, idx, road_wind_dt,
                                           beta_wind = 1, use_u_ref = TRUE, eps = 1e-12) {
  # attach grid_id to pm values
  pm_dt2 <- merge(pm_dt[, .(lon, lat, pm12 = val)], pmgrid, by = c("lon","lat"), all.x = TRUE)
  if (any(is.na(pm_dt2$grid_id))) stop("Some CMAQ cells could not be matched to pmgrid by lon/lat.")
  setkey(pm_dt2, grid_id)

  # road -> grid_id via idx
  dt <- copy(road_wind_dt)
  dt[, grid_id := pmgrid$grid_id[idx]]
  setkey(dt, grid_id)

  # attach pm12 by grid_id
  dt <- pm_dt2[dt]

  # roadiness relative factor (within-cell)
  dt[, road_eff := pmax(road, 0)]
  dt[, road_mean := mean(road_eff, na.rm = TRUE), by = grid_id]
  dt[, road_rel  := (road_eff + eps) / (road_mean + eps)]

  # meteo factor (wind speed only)
  if (isTRUE(use_u_ref)) {
    u_ref <- as.numeric(stats::median(dt$wspd, na.rm = TRUE))
    if (!is.finite(u_ref) || u_ref <= 0) u_ref <- 1
    dt[, meteo := ((u_ref + eps) / (wspd + eps))^beta_wind]
  } else {
    dt[, meteo := (inv_u)^beta_wind]
  }

  # combined modifier S_i (AADT will be multiplied later)
  dt[, S := road_rel * meteo]

  # RAW
  dt[, pm1_raw := pm12 * S]

  # NORM (keep cell mean)
  dt[, S_bar := mean(S, na.rm = TRUE), by = grid_id]
  dt[, pm1_norm := pm12 * (S / (S_bar + eps))]

  dt[, .(x, y, pm1_raw, pm1_norm)]
}

save_fst_xyval <- function(dt, outfile) {
  write_fst(dt, outfile)
  cat("WROTE: ", outfile, " | rows=", nrow(dt), " | size=", file.info(outfile)$size, "\n", sep="")
}

read_stack_years_scn <- function(layer = c("NRD","ONR","TOTAL"), years, metric, scenario = c("RAW","NORM")) {
  layer <- match.arg(layer)
  scenario <- match.arg(scenario)

  files <- switch(layer,
    NRD   = file.path(out_data, sprintf("VA_1km_NRD_%s_%d_%s.fst", metric, years, scenario)),
    ONR   = file.path(out_data, sprintf("VA_1km_ONR_%s_%d_%s.fst", metric, years, scenario)),
    TOTAL = file.path(out_data, sprintf("VA_1km_TOTAL_%s_%d_%s.fst", metric, years, scenario))
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

plot_panel_va <- function(stack_dt, out_pdf) {
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
# RUN
# =========================================================
cat("=== VA 1km downscale: roadiness + meteo (wind speed) | RAW + NORM ===\n")
cat("Output root: ", out_dir, "\n", sep="")

# load roadiness
road_nrd <- load_roadiness(road_col_NRD)
road_onr <- load_roadiness(road_col_ONR)

# build sf points once for NN mapping
road_sf_ref <- st_as_sf(road_nrd[, .(x,y)], coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

for (yr in years) {
  cat("\n-----------------------------\n")
  cat("YEAR: ", yr, " | metric: ", metric, "\n", sep="")
  cat("-----------------------------\n")

  # Read annual CMAQ
  nrd12 <- read_annual_rds("NRD", metric, yr)
  onr12 <- read_annual_rds("ONR", metric, yr)

  # Build pm grid + NN mapping index
  mm <- make_pmgrid_and_mapidx(nrd12, road_sf_ref)
  pmgrid <- mm$pmgrid
  idx    <- mm$idx

  # Extract wind to road grid
  road_nrd_w <- extract_wind_to_road(road_nrd, yr, p4s, eps)
  road_onr_w <- extract_wind_to_road(road_onr, yr, p4s, eps)

  # Downscale NRD
  nrd1 <- downscale_one_layer_road_meteo(
    nrd12, pmgrid, idx, road_nrd_w,
    beta_wind = beta_wind, use_u_ref = use_u_ref, eps = eps
  )

  f_nrd_raw  <- file.path(out_data, sprintf("VA_1km_NRD_%s_%d_RAW.fst",  metric, yr))
  f_nrd_norm <- file.path(out_data, sprintf("VA_1km_NRD_%s_%d_NORM.fst", metric, yr))
  save_fst_xyval(nrd1[, .(x, y, pm1 = pm1_raw)],  f_nrd_raw)
  save_fst_xyval(nrd1[, .(x, y, pm1 = pm1_norm)], f_nrd_norm)

  # Downscale ONR
  onr1 <- downscale_one_layer_road_meteo(
    onr12, pmgrid, idx, road_onr_w,
    beta_wind = beta_wind, use_u_ref = use_u_ref, eps = eps
  )

  f_onr_raw  <- file.path(out_data, sprintf("VA_1km_ONR_%s_%d_RAW.fst",  metric, yr))
  f_onr_norm <- file.path(out_data, sprintf("VA_1km_ONR_%s_%d_NORM.fst", metric, yr))
  save_fst_xyval(onr1[, .(x, y, pm1 = pm1_raw)],  f_onr_raw)
  save_fst_xyval(onr1[, .(x, y, pm1 = pm1_norm)], f_onr_norm)

  # TOTAL = NRD + ONR (RAW)
  nrd_raw <- nrd1[, .(x, y, nrd = pm1_raw)]
  onr_raw <- onr1[, .(x, y, onr = pm1_raw)]
  setkey(nrd_raw, x, y)
  setkey(onr_raw, x, y)
  tot_raw <- nrd_raw[onr_raw]
  tot_raw[, pm1_total := nrd + onr]
  tot_raw <- tot_raw[, .(x, y, pm1_total)]

  f_tot_raw <- file.path(out_data, sprintf("VA_1km_TOTAL_%s_%d_RAW.fst", metric, yr))
  save_fst_xyval(tot_raw, f_tot_raw)

  # TOTAL = NRD + ONR (NORM)
  nrd_norm <- nrd1[, .(x, y, nrd = pm1_norm)]
  onr_norm <- onr1[, .(x, y, onr = pm1_norm)]
  setkey(nrd_norm, x, y)
  setkey(onr_norm, x, y)
  tot_norm <- nrd_norm[onr_norm]
  tot_norm[, pm1_total := nrd + onr]
  tot_norm <- tot_norm[, .(x, y, pm1_total)]

  f_tot_norm <- file.path(out_data, sprintf("VA_1km_TOTAL_%s_%d_NORM.fst", metric, yr))
  save_fst_xyval(tot_norm, f_tot_norm)

  # quick diagnostics
  cat("WindSpeed_eff range: ", paste(range(road_nrd_w$wspd, na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("NRD RAW range:  ", paste(range(nrd1$pm1_raw,  na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("NRD NORM range: ", paste(range(nrd1$pm1_norm, na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("ONR RAW range:  ", paste(range(onr1$pm1_raw,  na.rm=TRUE), collapse=" ~ "), "\n", sep="")
  cat("ONR NORM range: ", paste(range(onr1$pm1_norm, na.rm=TRUE), collapse=" ~ "), "\n", sep="")

  rm(nrd12, onr12, mm, pmgrid, idx, road_nrd_w, road_onr_w,
     nrd1, onr1, nrd_raw, onr_raw, tot_raw, nrd_norm, onr_norm, tot_norm)
  gc()
}

cat("\n=== MAKING PANEL MAPS (2011–2020) | RAW vs NORM ===\n")

for (scn in c("RAW","NORM")) {
  stack_nrd <- read_stack_years_scn("NRD", years, metric, scenario = scn)
  pdf_nrd <- file.path(out_fig, sprintf("PANEL_VA_NRD_%s_2011_2020_%s.pdf", metric, scn))
  plot_panel_va(stack_nrd, pdf_nrd)

  stack_onr <- read_stack_years_scn("ONR", years, metric, scenario = scn)
  pdf_onr <- file.path(out_fig, sprintf("PANEL_VA_ONR_%s_2011_2020_%s.pdf", metric, scn))
  plot_panel_va(stack_onr, pdf_onr)

  stack_tot <- read_stack_years_scn("TOTAL", years, metric, scenario = scn)
  pdf_tot <- file.path(out_fig, sprintf("PANEL_VA_TOTAL_%s_2011_2020_%s.pdf", metric, scn))
  plot_panel_va(stack_tot, pdf_tot)

  rm(stack_nrd, stack_onr, stack_tot); gc()
}

cat("\nALL DONE ✅\n")
cat("Final data (FST): ", out_data, "\n", sep="")
cat("Panel maps (PDF): ", out_fig,  "\n", sep="")
