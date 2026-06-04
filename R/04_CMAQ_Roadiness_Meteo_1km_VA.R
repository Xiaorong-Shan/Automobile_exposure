#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(fst)
  library(raster)
  library(ggplot2)
  library(USAboundaries)
  library(viridis)
  library(scales)
})

# =========================================================
# SETTINGS
# =========================================================
state_abbr <- "VA"
years <- 2011:2020
metric <- "mean"

p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

base_dir   <- "/scratch/xshan2/R_Code/Automobiles"
annual_dir <- file.path(base_dir, "FIGURES")

roadiness_fst <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/VA/roadiness_1km_hw_loc_VA.fst"

road_col_NRD <- "nroad.dist2_loc"  # gasoline / local roads
road_col_ONR <- "nroad.dist2_hw"   # diesel / highway roads

wind_dir <- "/scratch/xshan2/NLDAS_wind/VA"

out_dir      <- file.path(base_dir, "DOWNSCALE_1KM_VA_WIND_2011_2020_SAVE_ALL")
out_weight   <- file.path(out_dir, "01_WIND_WEIGHTS_FST")
out_pm       <- file.path(out_dir, "02_PM1_WIND_FST")
out_total    <- file.path(out_dir, "03_TOTAL_FST")
out_fig      <- file.path(out_dir, "04_PANEL_FIGURES")
out_diag     <- file.path(out_dir, "05_DIAGNOSTICS")

dir.create(out_weight, recursive = TRUE, showWarnings = FALSE)
dir.create(out_pm,     recursive = TRUE, showWarnings = FALSE)
dir.create(out_total,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig,    recursive = TRUE, showWarnings = FALSE)
dir.create(out_diag,   recursive = TRUE, showWarnings = FALSE)

eps <- 1e-12
U0 <- 0.1
beta_wind <- 0.5

legend_upper_q <- 0.99
facet_ncol <- 5

# =========================================================
# HELPERS
# =========================================================
msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

get_va_boundary <- function() {
  va <- USAboundaries::us_states(resolution = "low")
  va <- va[va$state_abbr == state_abbr, ]
  st_transform(st_as_sf(va), st_crs(p4s))
}

read_annual_rds <- function(sector, metric, yr) {
  f <- file.path(annual_dir, sprintf("ANNUAL_%s_%s_%d.rds", sector, metric, yr))
  if (!file.exists(f)) stop("Missing annual RDS: ", f)

  dt <- readRDS(f)
  setDT(dt)

  if (!all(c("lon", "lat", "val") %in% names(dt))) {
    stop("Annual RDS must contain lon, lat, val: ", f)
  }

  dt <- dt[is.finite(val) & !is.na(val)]
  dt
}

load_roadiness <- function(colname) {
  rd <- read_fst(roadiness_fst, as.data.table = TRUE)
  setDT(rd)

  if (!all(c("x", "y") %in% names(rd))) stop("Roadiness fst must contain x,y.")
  if (!(colname %in% names(rd))) stop("Missing roadiness column: ", colname)

  dt <- rd[, .(
    x = x,
    y = y,
    road = as.numeric(get(colname))
  )]

  dt[!is.finite(road) | is.na(road), road := 0]
  dt[road < 0, road := 0]

  dt
}

detect_grid_spacing <- function(dt) {
  list(
    dx = median(diff(sort(unique(dt$x))), na.rm = TRUE),
    dy = median(diff(sort(unique(dt$y))), na.rm = TRUE)
  )
}

extract_wind <- function(road_dt, yr) {
  msg("Extracting wind for year:", yr)

  wind_tif <- file.path(wind_dir, sprintf("NLDAS_%s_%d_Wind_EN_mean.tif", state_abbr, yr))
  if (!file.exists(wind_tif)) stop("Missing wind tif: ", wind_tif)

  wd <- raster::brick(wind_tif)

  pts <- st_as_sf(
    road_dt[, .(x, y)],
    coords = c("x", "y"),
    crs = st_crs(p4s),
    remove = FALSE
  )

  pts_ll <- st_transform(pts, crs(wd))
  xy <- st_coordinates(pts_ll)

  u <- raster::extract(wd[[1]], xy)
  v <- raster::extract(wd[[2]], xy)

  dt <- copy(road_dt)
  dt[, u := u]
  dt[, v := v]
  dt[, windspeed := sqrt(u^2 + v^2)]
  dt[, theta_to := atan2(v, u)]

  ws_med <- median(dt$windspeed, na.rm = TRUE)
  dt[!is.finite(windspeed) | is.na(windspeed), windspeed := ws_med]
  dt[!is.finite(theta_to) | is.na(theta_to), theta_to := 0]

  dt
}

attach_neighbor_road <- function(dt, name, sx, sy, dx_grid, dy_grid) {
  tmp <- copy(dt)

  tmp[, nx := x + sx * dx_grid]
  tmp[, ny := y + sy * dy_grid]

  src <- dt[, .(
    nx = x,
    ny = y,
    road_neighbor = road
  )]

  setkey(src, nx, ny)
  setkey(tmp, nx, ny)

  tmp <- src[tmp]

  new_col <- paste0("road_", name)
  setnames(tmp, "road_neighbor", new_col)

  tmp[!is.finite(get(new_col)) | is.na(get(new_col)), (new_col) := road]
  tmp[, c("nx", "ny") := NULL]

  tmp
}

make_wind_weight <- function(road_wind_dt) {
  msg("Computing 8-neighbor upwind roadiness + wind-speed beta factor ...")

  dt <- copy(road_wind_dt)

  gs <- detect_grid_spacing(dt)
  dx_grid <- gs$dx
  dy_grid <- gs$dy

  msg("Grid spacing:", dx_grid, dy_grid)

  dt <- attach_neighbor_road(dt, "E",   1,  0, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "NE",  1,  1, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "N",   0,  1, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "NW", -1,  1, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "W",  -1,  0, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "SW", -1, -1, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "S",   0, -1, dx_grid, dy_grid)
  dt <- attach_neighbor_road(dt, "SE",  1, -1, dx_grid, dy_grid)

  # Upwind source direction = opposite of wind blowing-to direction
  dt[, theta_upwind := (theta_to + pi) %% (2 * pi)]

  dirs <- data.table(
    name  = c("E", "NE", "N", "NW", "W", "SW", "S", "SE"),
    angle = c(0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4),
    col   = c("road_E", "road_NE", "road_N", "road_NW", "road_W", "road_SW", "road_S", "road_SE")
  )

  for (k in seq_len(nrow(dirs))) {
    cname <- paste0("angdiff_", dirs$name[k])
    wname <- paste0("w_", dirs$name[k])
    a <- dirs$angle[k]

    dt[, (cname) := abs(atan2(sin(theta_upwind - a), cos(theta_upwind - a)))]
    dt[, (wname) := 1 / (get(cname) + 0.05)]
  }

  weight_cols <- paste0("w_", dirs$name)
  road_cols   <- dirs$col

  dt[, w_sum := rowSums(.SD, na.rm = TRUE), .SDcols = weight_cols]
  dt[, road_upwind8 := 0]

  for (k in seq_along(weight_cols)) {
    dt[, road_upwind8 := road_upwind8 + get(weight_cols[k]) * get(road_cols[k])]
  }

  dt[, road_upwind8 := road_upwind8 / (w_sum + eps)]
  dt[!is.finite(road_upwind8) | is.na(road_upwind8), road_upwind8 := road]

  ws_mean <- mean(dt$windspeed, na.rm = TRUE)

  dt[, speed_factor_beta := ((ws_mean + eps) / (windspeed + U0))^beta_wind]

  # final wind-adjusted redistribution weight
  dt[, W_wind := road_upwind8 * speed_factor_beta]

  dt[!is.finite(W_wind) | is.na(W_wind), W_wind := 0]
  dt[W_wind < 0, W_wind := 0]

  # keep useful intermediate variables
  out <- dt[, .(
    x, y,
    road,
    u, v,
    windspeed,
    theta_to,
    theta_upwind,
    road_upwind8,
    speed_factor_beta,
    W_wind
  )]

  out
}

make_pmgrid_and_mapidx <- function(cmaq_dt, road_dt) {
  pmgrid <- unique(cmaq_dt[, .(lon, lat)])
  pmgrid[, grid_id := .I]

  pm_sf <- st_as_sf(
    pmgrid,
    coords = c("lon", "lat"),
    crs = 4326,
    remove = FALSE
  )
  pm_sf <- st_transform(pm_sf, st_crs(p4s))

  road_sf <- st_as_sf(
    road_dt[, .(x, y)],
    coords = c("x", "y"),
    crs = st_crs(p4s),
    remove = FALSE
  )

  idx <- st_nearest_feature(road_sf, pm_sf)

  list(pmgrid = pmgrid, idx = idx)
}

downscale_with_weight <- function(cmaq_dt, pmgrid, idx, weight_dt) {
  pm_dt2 <- merge(
    cmaq_dt[, .(lon, lat, pm12 = val)],
    pmgrid,
    by = c("lon", "lat"),
    all.x = TRUE
  )

  if (any(is.na(pm_dt2$grid_id))) {
    stop("Some CMAQ cells could not be matched to pmgrid.")
  }

  setkey(pm_dt2, grid_id)

  dt <- copy(weight_dt)
  dt[, grid_id := pmgrid$grid_id[idx]]

  setkey(dt, grid_id)
  dt <- pm_dt2[dt]

  dt[!is.finite(W_wind) | is.na(W_wind), W_wind := 0]
  dt[W_wind < 0, W_wind := 0]

  # Mass conservation within each CMAQ grid:
  # mean(pm1 within grid_id) = pm12
  dt[, W_mean_12km := mean(W_wind, na.rm = TRUE), by = grid_id]
  dt[, weight_norm := (W_wind + eps) / (W_mean_12km + eps)]
  dt[, pm1 := pm12 * weight_norm]

  check <- dt[, .(
    pm12 = unique(pm12)[1],
    pm1_mean = mean(pm1, na.rm = TRUE),
    diff = mean(pm1, na.rm = TRUE) - unique(pm12)[1]
  ), by = grid_id]

  max_error <- max(abs(check$diff), na.rm = TRUE)
  msg("Max conservation error:", signif(max_error, 6))

  pm_out <- dt[, .(
    x, y,
    grid_id,
    pm12,
    road,
    road_upwind8,
    windspeed,
    speed_factor_beta,
    W_wind,
    W_mean_12km,
    weight_norm,
    pm1
  )]

  diag_out <- data.table(
    n_grid = nrow(check),
    max_abs_conservation_error = max_error,
    pm1_min = min(pm_out$pm1, na.rm = TRUE),
    pm1_mean = mean(pm_out$pm1, na.rm = TRUE),
    pm1_max = max(pm_out$pm1, na.rm = TRUE),
    W_min = min(pm_out$W_wind, na.rm = TRUE),
    W_mean = mean(pm_out$W_wind, na.rm = TRUE),
    W_max = max(pm_out$W_wind, na.rm = TRUE)
  )

  list(pm = pm_out, diag = diag_out)
}

save_fst_log <- function(dt, f) {
  write_fst(dt, f)
  msg("WROTE:", f, "| rows =", nrow(dt), "| size =", file.info(f)$size)
}

read_stack_years <- function(layer, years) {
  files <- switch(
    layer,
    NRD = file.path(out_pm, sprintf("VA_1km_NRD_wind_%s_%d.fst", metric, years)),
    ONR = file.path(out_pm, sprintf("VA_1km_ONR_wind_%s_%d.fst", metric, years)),
    TOTAL = file.path(out_total, sprintf("VA_1km_TOTAL_wind_%s_%d.fst", metric, years))
  )

  lst <- lapply(seq_along(years), function(i) {
    f <- files[i]
    if (!file.exists(f)) stop("Missing file for panel: ", f)

    dt <- read_fst(f, as.data.table = TRUE)
    setDT(dt)

    if (layer == "TOTAL") {
      setnames(dt, "pm1_total", "val")
    } else {
      setnames(dt, "pm1", "val")
    }

    dt[, year := years[i]]
    dt[, .(x, y, val, year)]
  })

  rbindlist(lst, use.names = TRUE)
}

plot_panel_va <- function(stack_dt, out_pdf, title_text) {
  va_sf <- get_va_boundary()

  stack_dt <- stack_dt[is.finite(val) & !is.na(val)]
  stack_dt[, year_f := factor(year, levels = sort(unique(year)))]

  xlim <- range(stack_dt$x, na.rm = TRUE)
  ylim <- range(stack_dt$y, na.rm = TRUE)

  lim_upper <- as.numeric(quantile(stack_dt$val, legend_upper_q, na.rm = TRUE))
  if (!is.finite(lim_upper) || lim_upper <= 0) {
    lim_upper <- max(stack_dt$val, na.rm = TRUE)
  }

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
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")
    )

  grDevices::pdf(out_pdf, width = 16, height = 6.5, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()

  msg("WROTE:", out_pdf, "| size =", file.info(out_pdf)$size)
}

# =========================================================
# RUN
# =========================================================
msg("=== VA 2011-2020 wind-adjusted 1km downscaling ===")
msg("beta_wind =", beta_wind, "| U0 =", U0)
msg("Output root:", out_dir)

road_nrd_base <- load_roadiness(road_col_NRD)
road_onr_base <- load_roadiness(road_col_ONR)

all_diag <- list()

for (yr in years) {
  msg("--------------------------------------------------")
  msg("YEAR:", yr)
  msg("--------------------------------------------------")

  # -------------------------
  # Read CMAQ annual fields
  # -------------------------
  nrd12 <- read_annual_rds("NRD", metric, yr)
  onr12 <- read_annual_rds("ONR", metric, yr)

  # -------------------------
  # Wind weights: NRD
  # -------------------------
  msg("NRD / gasoline wind weight ...")
  nrd_wind_input <- extract_wind(road_nrd_base, yr)
  nrd_weight <- make_wind_weight(nrd_wind_input)
  nrd_weight[, sector := "NRD"]
  nrd_weight[, year := yr]

  f_nrd_w <- file.path(out_weight, sprintf("VA_1km_NRD_wind_weight_%s_%d.fst", metric, yr))
  save_fst_log(nrd_weight, f_nrd_w)

  # -------------------------
  # Wind weights: ONR
  # -------------------------
  msg("ONR / diesel wind weight ...")
  onr_wind_input <- extract_wind(road_onr_base, yr)
  onr_weight <- make_wind_weight(onr_wind_input)
  onr_weight[, sector := "ONR"]
  onr_weight[, year := yr]

  f_onr_w <- file.path(out_weight, sprintf("VA_1km_ONR_wind_weight_%s_%d.fst", metric, yr))
  save_fst_log(onr_weight, f_onr_w)

  # -------------------------
  # CMAQ mapping
  # Use NRD grid as reference; ONR should share same grid
  # -------------------------
  msg("Building CMAQ mapping ...")
  mm <- make_pmgrid_and_mapidx(nrd12, road_nrd_base)
  pmgrid <- mm$pmgrid
  idx <- mm$idx

  # -------------------------
  # Downscale NRD
  # -------------------------
  msg("Downscaling NRD / gasoline ...")
  nrd_res <- downscale_with_weight(nrd12, pmgrid, idx, nrd_weight)
  nrd1 <- nrd_res$pm
  nrd1[, sector := "NRD"]
  nrd1[, year := yr]

  f_nrd_pm <- file.path(out_pm, sprintf("VA_1km_NRD_wind_%s_%d.fst", metric, yr))
  save_fst_log(nrd1, f_nrd_pm)

  nrd_diag <- copy(nrd_res$diag)
  nrd_diag[, sector := "NRD"]
  nrd_diag[, year := yr]

  # -------------------------
  # Downscale ONR
  # -------------------------
  msg("Downscaling ONR / diesel ...")
  onr_res <- downscale_with_weight(onr12, pmgrid, idx, onr_weight)
  onr1 <- onr_res$pm
  onr1[, sector := "ONR"]
  onr1[, year := yr]

  f_onr_pm <- file.path(out_pm, sprintf("VA_1km_ONR_wind_%s_%d.fst", metric, yr))
  save_fst_log(onr1, f_onr_pm)

  onr_diag <- copy(onr_res$diag)
  onr_diag[, sector := "ONR"]
  onr_diag[, year := yr]

  # -------------------------
  # TOTAL = NRD + ONR
  # -------------------------
  msg("Creating TOTAL = NRD + ONR ...")

  nrd_tmp <- nrd1[, .(x, y, pm1_NRD = pm1)]
  onr_tmp <- onr1[, .(x, y, pm1_ONR = pm1)]

  setkey(nrd_tmp, x, y)
  setkey(onr_tmp, x, y)

  total1 <- nrd_tmp[onr_tmp]
  total1[, pm1_total := pm1_NRD + pm1_ONR]
  total1[, year := yr]
  total1 <- total1[, .(x, y, pm1_NRD, pm1_ONR, pm1_total, year)]

  f_total <- file.path(out_total, sprintf("VA_1km_TOTAL_wind_%s_%d.fst", metric, yr))
  save_fst_log(total1, f_total)

  total_diag <- data.table(
    sector = "TOTAL",
    year = yr,
    n_grid = NA_integer_,
    max_abs_conservation_error = NA_real_,
    pm1_min = min(total1$pm1_total, na.rm = TRUE),
    pm1_mean = mean(total1$pm1_total, na.rm = TRUE),
    pm1_max = max(total1$pm1_total, na.rm = TRUE),
    W_min = NA_real_,
    W_mean = NA_real_,
    W_max = NA_real_
  )

  all_diag[[paste0("NRD_", yr)]] <- nrd_diag
  all_diag[[paste0("ONR_", yr)]] <- onr_diag
  all_diag[[paste0("TOTAL_", yr)]] <- total_diag

  rm(nrd12, onr12, nrd_wind_input, onr_wind_input,
     nrd_weight, onr_weight, nrd_res, onr_res,
     nrd1, onr1, total1, mm, pmgrid, idx)
  gc()
}

# =========================================================
# SAVE DIAGNOSTICS
# =========================================================
diag_dt <- rbindlist(all_diag, use.names = TRUE, fill = TRUE)
diag_csv <- file.path(out_diag, "DIAGNOSTICS_wind_downscale_2011_2020.csv")
fwrite(diag_dt, diag_csv)
msg("WROTE:", diag_csv)

# =========================================================
# PANEL MAPS
# =========================================================
msg("Making panel maps ...")

stack_nrd <- read_stack_years("NRD", years)
stack_onr <- read_stack_years("ONR", years)
stack_tot <- read_stack_years("TOTAL", years)

plot_panel_va(
  stack_nrd,
  file.path(out_fig, sprintf("PANEL_VA_NRD_gasoline_wind_%s_2011_2020.pdf", metric)),
  "VA 1 km NRD / gasoline wind-adjusted PM2.5, 2011-2020"
)

plot_panel_va(
  stack_onr,
  file.path(out_fig, sprintf("PANEL_VA_ONR_diesel_wind_%s_2011_2020.pdf", metric)),
  "VA 1 km ONR / diesel wind-adjusted PM2.5, 2011-2020"
)

plot_panel_va(
  stack_tot,
  file.path(out_fig, sprintf("PANEL_VA_TOTAL_wind_%s_2011_2020.pdf", metric)),
  "VA 1 km TOTAL wind-adjusted PM2.5, 2011-2020"
)

msg("ALL DONE")
msg("Output root:", out_dir)
