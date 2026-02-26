#!/usr/bin/env Rscript
# ============================================================
# Mean-conserving 1 km downscaling of 12 km CMAQ annual PM2.5
# using:
#   C_i = C_g * [ r_i * U_i^{-1} * f(Δθ_i) ] / mean_{j in g}[ r_j * U_j^{-1} * f(Δθ_j) ]
#
# Key fixes in this version:
#   (A) CLIP the 1 km grid to the state boundary (so "VA" is truly VA-only).
#       This prevents massive NA wind extraction when wind tif is state-masked.
#   (B) Wind extraction transforms points to the wind raster CRS (no hard-coded lon/lat).
#   (C) Panel maps are skipped if yearly FSTs are missing.
# ============================================================

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

# ============================================================
# SETTINGS (EDIT FOR EACH STATE)
# ============================================================

# ---- State ID
state_abbr <- "VA"

# ---- Projection used by your 1 km roadiness grid (LCC)
p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

# ---- Base directories
base_dir   <- "/scratch/xshan2/R_Code/Automobiles"
annual_dir <- file.path(base_dir, "FIGURES")  # contains ANNUAL_{NRD,ONR}_{mean|median}_{year}.rds

# ---- 1 km roadiness grid (must contain x,y and roadiness columns)
roadiness_fst <- sprintf("/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/%s/roadiness_1km_hw_loc_%s.fst",
                         state_abbr, state_abbr)

# ---- Road centerlines (for computing φ_i once)
roads_shp <- "/scratch/xshan2/R_Code/Automobiles/roadiness/TRAN_Virginia_State_Shape/Shape/Trans_RoadSegment.shp"

# Optional: filter expression to keep only major roads (set NULL to keep all).
roads_filter_expr <- NULL

# ---- Wind GeoTIFF: annual mean u/v (Band1=Wind_E=u, Band2=Wind_N=v)
wind_dir     <- sprintf("/scratch/xshan2/NLDAS_wind/%s", state_abbr)
wind_tif_fmt <- file.path(wind_dir, sprintf("NLDAS_%s_%%d_Wind_EN_mean.tif", state_abbr))

# ---- Years and CMAQ metric
years  <- 2011:2020
metric <- "mean"  # or "median"

# ---- Roadiness columns
road_col_NRD <- "nroad.dist2_loc"
road_col_ONR <- "nroad.dist2_hw"

# ---- Alignment function choice
# "cos01"  : f(Δθ) = (1 + cos(Δθ))/2   in [0,1]
# "cospos" : f(Δθ) = max(cos(Δθ), 0)   in [0,1]
align_fun <- "cos01"

# ---- Chunk size for nearest-road matching (φ computation)
phi_chunk_size <- 50000

# ---- Numerical safeguard for division by zero only
eps <- 1e-12

# ---- Output folders
out_dir  <- file.path(base_dir, sprintf("DOWNSCALE_1KM_%s_ROAD_WINDDIR_NORMONLY", state_abbr))
out_data <- file.path(out_dir, "DATA_FST")
out_fig  <- file.path(out_dir, "FIGURES")
out_aux  <- file.path(out_dir, "AUX")  # cache phi here
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_aux,  recursive = TRUE, showWarnings = FALSE)

phi_cache_fst <- file.path(out_aux, sprintf("%s_1km_phi_from_roads.fst", state_abbr))

# ---- Panel plot options (optional)
make_panel_maps <- TRUE
facet_ncol <- 5
legend_upper_q <- 0.99

# ============================================================
# HELPERS
# ============================================================

msg <- function(...) cat(paste0(..., "\n"))

stopif_has_cols <- function(dt, cols, where = "object") {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) stop(where, " is missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

# ------------------------------------------------------------
# Clip 1 km points to state boundary (keep only points within state)
# This is the MOST important fix if your 1 km grid is a rectangle bbox
# and your wind tif is masked to the state boundary.
# ------------------------------------------------------------
clip_points_to_state <- function(dt_xy, state_abbr, p4s) {
  stopif_has_cols(dt_xy, c("x","y"), where = "dt_xy for clip")

  st_outline <- USAboundaries::us_states(resolution = "low")
  st_outline <- st_outline[st_outline$state_abbr == state_abbr, ]
  if (nrow(st_outline) == 0) stop("State not found in USAboundaries: ", state_abbr)

  st_outline <- st_transform(st_as_sf(st_outline), st_crs(p4s))

  pts <- st_as_sf(dt_xy, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)
  inside <- st_within(pts, st_outline, sparse = FALSE)[,1]

  msg("Clip to ", state_abbr, ": kept ", sum(inside), " / ", nrow(dt_xy),
      " (", round(100*mean(inside), 2), "%) points")

  dt_xy[inside]
}

# ------------------------------------------------------------
# Read annual CMAQ RDS: lon, lat, val
# ------------------------------------------------------------
read_annual_rds <- function(sector, metric, yr) {
  f <- file.path(annual_dir, sprintf("ANNUAL_%s_%s_%d.rds", sector, metric, yr))
  if (!file.exists(f)) stop("Missing annual RDS: ", f)

  dt <- readRDS(f)
  setDT(dt)

  stopif_has_cols(dt, c("lon","lat","val"), where = paste0("Annual RDS ", f))
  dt <- dt[is.finite(val) & !is.na(val)]
  dt
}

# ------------------------------------------------------------
# Load 1 km roadiness grid with a selected column -> (x,y,road)
# ------------------------------------------------------------
load_roadiness <- function(colname) {
  if (!file.exists(roadiness_fst)) stop("roadiness_fst not found: ", roadiness_fst)

  rd <- read_fst(roadiness_fst, as.data.table = TRUE)
  setDT(rd)

  stopif_has_cols(rd, c("x","y"), where = "roadiness_fst")
  if (!(colname %in% names(rd))) stop("Roadiness col not found: ", colname)

  out <- rd[, .(x, y, road = as.numeric(get(colname)))]
  out[!is.finite(road) | is.na(road), road := 0]
  out
}

# ------------------------------------------------------------
# Build unique CMAQ centroid grid and map each 1 km cell to its
# nearest CMAQ centroid (index vector length = n_1km)
# ------------------------------------------------------------
make_pmgrid_and_mapidx <- function(pm_dt, road_sf_ref, p4s) {
  pmgrid <- unique(pm_dt[, .(lon, lat)])
  pmgrid[, grid_id := .I]

  pm_sf <- st_as_sf(pmgrid, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  pm_sf <- st_transform(pm_sf, st_crs(p4s))

  idx <- st_nearest_feature(road_sf_ref, pm_sf)
  list(pmgrid = pmgrid, idx = idx)
}

# ------------------------------------------------------------
# Compute per-road bearing (phi) in radians:
# 0=East, pi/2=North, range [-pi, pi]
# ------------------------------------------------------------
bearing_from_linestring <- function(geom) {
  coords <- st_coordinates(geom)
  if (nrow(coords) < 2) return(NA_real_)
  x1 <- coords[1, "X"];  y1 <- coords[1, "Y"]
  x2 <- coords[nrow(coords), "X"]; y2 <- coords[nrow(coords), "Y"]
  atan2(y2 - y1, x2 - x1)
}

# ------------------------------------------------------------
# Compute φ_i for each 1 km grid cell i:
# φ_i = bearing of nearest road segment
# Cached to phi_cache_fst for reuse across years.
# ------------------------------------------------------------
compute_or_load_phi <- function(road_xy_dt,
                                roads_shp,
                                p4s,
                                phi_cache_fst,
                                chunk_size = 50000,
                                filter_expr = NULL) {
  stopif_has_cols(road_xy_dt, c("x","y"), where = "road_xy_dt for phi")

  if (file.exists(phi_cache_fst)) {
    msg("Loading cached phi: ", phi_cache_fst)
    phi_dt <- read_fst(phi_cache_fst, as.data.table = TRUE)
    setDT(phi_dt)
    stopif_has_cols(phi_dt, c("x","y","phi"), where = "phi_cache_fst")
    return(phi_dt)
  }

  if (!file.exists(roads_shp)) stop("roads_shp not found: ", roads_shp)

  msg("Reading roads: ", roads_shp)
  roads <- st_read(roads_shp, quiet = TRUE)
  if (nrow(roads) == 0) stop("roads_shp has 0 features: ", roads_shp)

  if (!is.null(filter_expr)) {
    roads <- roads[eval(filter_expr, envir = roads), ]
    if (nrow(roads) == 0) stop("No roads left after filtering. Check roads_filter_expr.")
  }

  roads <- st_transform(roads, st_crs(p4s))
  roads <- st_cast(roads, "LINESTRING", warn = FALSE)

  msg("Computing road bearings (phi) for ", nrow(roads), " segments ...")
  roads$phi <- vapply(st_geometry(roads), bearing_from_linestring, numeric(1))
  if (sum(!is.finite(roads$phi)) > 0) {
    stop("Some road bearings are NA/Inf. Check road geometry.")
  }

  pts <- st_as_sf(road_xy_dt, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

  n <- nrow(pts)
  idx <- integer(n)

  msg("Computing nearest-road index for ", n, " grid points (chunk_size=", chunk_size, ") ...")
  t0 <- Sys.time()
  for (s in seq(1, n, by = chunk_size)) {
    e <- min(s + chunk_size - 1, n)
    idx[s:e] <- st_nearest_feature(pts[s:e, ], roads)
    msg("  done: ", e, " / ", n)
  }
  msg("Elapsed (min): ", round(difftime(Sys.time(), t0, units = "mins"), 2))

  phi_dt <- copy(road_xy_dt)
  phi_dt[, phi := as.numeric(roads$phi[idx])]

  if (sum(!is.finite(phi_dt$phi)) > 0) stop("phi_dt contains NA/Inf, unexpected.")

  write_fst(phi_dt[, .(x,y,phi)], phi_cache_fst)
  msg("WROTE phi cache: ", phi_cache_fst, " | rows=", nrow(phi_dt))

  phi_dt[, .(x,y,phi)]
}

# ------------------------------------------------------------
# Extract annual mean wind (u/v) to 1 km centroids:
# - Band 1 = u (Wind_E)
# - Band 2 = v (Wind_N)
# Returns U and theta in radians
#
# IMPORTANT:
#   Transform points to the WIND RASTER CRS (not hard-coded lon/lat).
# ------------------------------------------------------------
extract_wind_uv_to_road <- function(road_xy_dt, yr, p4s, eps = 1e-12, method = "bilinear") {
  wind_tif <- sprintf(wind_tif_fmt, yr)
  if (!file.exists(wind_tif)) stop("Missing wind tif: ", wind_tif)

  wd <- raster::brick(wind_tif)
  if (raster::nlayers(wd) < 2) stop("Wind tif must have >=2 bands (u,v): ", wind_tif)

  rcrs_obj <- raster::crs(wd)
  rcrs_txt <- as.character(rcrs_obj)
  if (is.na(rcrs_txt) || nchar(rcrs_txt) == 0) stop("Wind raster has missing CRS: ", wind_tif)

  pts_lcc <- st_as_sf(road_xy_dt, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)
  pts_ras <- st_transform(pts_lcc, st_crs(rcrs_txt))
  xy      <- st_coordinates(pts_ras)

  u <- raster::extract(wd[[1]], xy, method = method)
  v <- raster::extract(wd[[2]], xy, method = method)

  out <- copy(road_xy_dt)
  out[, `:=`(u = as.numeric(u), v = as.numeric(v))]

  bad <- !is.finite(out$u) | !is.finite(out$v)
  if (any(bad)) {
    # After CLIPPING to state, this should be near zero.
    stop("Wind extraction still has non-finite u/v for ", sum(bad), " points. ",
         "This suggests the wind tif has NA even inside the state domain. ",
         "Check how the wind tif was created/masked.")
  }

  out[, `:=`(
    U     = sqrt(u^2 + v^2) + eps,
    theta = atan2(v, u)
  )]

  out
}

# ------------------------------------------------------------
# Alignment function f(Δθ)
# ------------------------------------------------------------
align_factor <- function(delta, mode = c("cos01","cospos")) {
  mode <- match.arg(mode)
  if (mode == "cos01")  return((1 + cos(delta)) / 2)
  if (mode == "cospos") return(pmax(cos(delta), 0))
}

# ------------------------------------------------------------
# Downscale one layer (NRD or ONR) using the formula above
# Returns: (x,y,pm1)
# ------------------------------------------------------------
downscale_one_layer_formula <- function(pm_dt, pmgrid, idx_cmaq, road_dt, eps = 1e-12, align_fun = "cos01") {

  stopif_has_cols(pm_dt, c("lon","lat","val"), where = "pm_dt")
  stopif_has_cols(pmgrid, c("lon","lat","grid_id"), where = "pmgrid")
  stopif_has_cols(road_dt, c("x","y","road","phi","U","theta"), where = "road_dt")

  pm_dt2 <- merge(pm_dt[, .(lon, lat, Cg = val)], pmgrid, by = c("lon","lat"), all.x = TRUE)
  if (any(is.na(pm_dt2$grid_id))) {
    stop("Some CMAQ cells could not be matched to pmgrid by lon/lat. Check grid consistency.")
  }
  setkey(pm_dt2, grid_id)

  dt <- copy(road_dt)
  dt[, grid_id := pmgrid$grid_id[idx_cmaq]]
  setkey(dt, grid_id)

  dt <- pm_dt2[dt]

  dt[, r := pmax(road, 0)]
  dt[, delta := atan2(sin(theta - phi), cos(theta - phi))]
  dt[, f := align_factor(delta, mode = align_fun)]
  dt[, S := r * (1 / U) * f]
  dt[, Sbar := mean(S, na.rm = TRUE), by = grid_id]
  dt[, Ci := Cg * (S / (Sbar + eps))]

  dt[, .(x, y, pm1 = Ci)]
}

# ------------------------------------------------------------
# Save FST
# ------------------------------------------------------------
save_fst <- function(dt, outfile) {
  write_fst(dt, outfile)
  msg("WROTE: ", outfile, " | rows=", nrow(dt), " | size=", file.info(outfile)$size)
}

# ------------------------------------------------------------
# Stack years for panel maps
# ------------------------------------------------------------
read_stack_years <- function(layer = c("NRD","ONR","TOTAL"), years, metric, out_data) {
  layer <- match.arg(layer)

  files <- switch(layer,
    NRD   = file.path(out_data, sprintf("%s_1km_NRD_%s_%d.fst",   state_abbr, metric, years)),
    ONR   = file.path(out_data, sprintf("%s_1km_ONR_%s_%d.fst",   state_abbr, metric, years)),
    TOTAL = file.path(out_data, sprintf("%s_1km_TOTAL_%s_%d.fst", state_abbr, metric, years))
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

# ------------------------------------------------------------
# Panel plot (optional)
# ------------------------------------------------------------
plot_panel_state <- function(stack_dt, out_pdf, state_abbr, p4s,
                             facet_ncol = 5, legend_upper_q = 0.99) {

  st_outline <- USAboundaries::us_states(resolution = "low")
  st_outline <- st_outline[st_outline$state_abbr == state_abbr, ]
  st_outline <- st_transform(st_as_sf(st_outline), st_crs(p4s))

  stack_dt <- stack_dt[is.finite(val) & !is.na(val)]
  stack_dt[, year_f := factor(year, levels = sort(unique(year)))]

  xlim <- range(stack_dt$x, na.rm = TRUE)
  ylim <- range(stack_dt$y, na.rm = TRUE)

  lim_upper <- as.numeric(quantile(stack_dt$val, legend_upper_q, na.rm = TRUE))
  if (!is.finite(lim_upper) || lim_upper <= 0) lim_upper <- max(stack_dt$val, na.rm = TRUE)

  p <- ggplot() +
    geom_raster(data = stack_dt, aes(x = x, y = y, fill = val)) +
    geom_sf(data = st_outline, fill = NA, color = "grey35", linewidth = 0.2) +
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
  msg("WROTE: ", out_pdf, " | size=", file.info(out_pdf)$size)
}

# ============================================================
# RUN
# ============================================================

msg("=== ", state_abbr, " | 1 km downscaling using r * U^{-1} * f(Δθ) | NORM ONLY ===")
msg("Output root: ", out_dir)

# 1) Load roadiness grids
road_nrd <- load_roadiness(road_col_NRD)
road_onr <- load_roadiness(road_col_ONR)

# 1.1) CRITICAL: clip to state boundary (make the grid truly state-only)
road_nrd <- clip_points_to_state(road_nrd, state_abbr, p4s)
road_onr <- clip_points_to_state(road_onr, state_abbr, p4s)

# 1.2) Align NRD/ONR to the exact same (x,y) set
setkey(road_nrd, x, y)
setkey(road_onr, x, y)
road_onr <- road_onr[road_nrd]   # keep only shared points
road_nrd <- road_nrd[road_onr]   # symmetric (safety)

if (nrow(road_nrd) == 0) stop("After clipping/alignment, 0 points remain. Check p4s and state boundary clip.")
msg("Final 1 km grid points: ", nrow(road_nrd))

# 2) Compute or load phi cache (use clipped grid!)
phi_dt <- compute_or_load_phi(
  road_xy_dt     = road_nrd[, .(x,y)],
  roads_shp      = roads_shp,
  p4s            = p4s,
  phi_cache_fst  = phi_cache_fst,
  chunk_size     = phi_chunk_size,
  filter_expr    = roads_filter_expr
)

# 3) Attach phi to both NRD and ONR roadiness tables
setkey(phi_dt, x, y)
setkey(road_nrd, x, y)
setkey(road_onr, x, y)

road_nrd <- phi_dt[road_nrd]
road_onr <- phi_dt[road_onr]

if (any(!is.finite(road_nrd$phi))) stop("phi missing in road_nrd after merge.")
if (any(!is.finite(road_onr$phi))) stop("phi missing in road_onr after merge.")

msg("phi attached. NA count (NRD): ", sum(!is.finite(road_nrd$phi)))

# 4) Build sf points once for CMAQ mapping
road_sf_ref <- st_as_sf(road_nrd[, .(x,y)], coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

# 5) Main loop over years
for (yr in years) {
  msg("\n-----------------------------")
  msg("YEAR: ", yr, " | metric: ", metric)
  msg("-----------------------------")

  # 5.1) Read annual CMAQ (12 km)
  nrd12 <- read_annual_rds("NRD", metric, yr)
  onr12 <- read_annual_rds("ONR", metric, yr)

  # 5.2) CMAQ mapping indices (based on NRD CMAQ grid)
  mm     <- make_pmgrid_and_mapidx(nrd12, road_sf_ref, p4s)
  pmgrid <- mm$pmgrid
  idx    <- mm$idx

  # 5.3) Wind extraction at 1 km grid (u/v -> U, theta)
  wind_nrd <- extract_wind_uv_to_road(road_nrd[, .(x,y)], yr, p4s, eps)
  wind_onr <- extract_wind_uv_to_road(road_onr[, .(x,y)], yr, p4s, eps)

  # 5.4) Build per-year road table inputs
  road_nrd_year <- road_nrd[, .(x,y,road,phi)]
  road_onr_year <- road_onr[, .(x,y,road,phi)]

  setkey(road_nrd_year, x, y)
  setkey(road_onr_year, x, y)
  setkey(wind_nrd, x, y)
  setkey(wind_onr, x, y)

  road_nrd_year <- wind_nrd[road_nrd_year]
  road_onr_year <- wind_onr[road_onr_year]

  if (any(!is.finite(road_nrd_year$U))) stop("Non-finite U in NRD wind extraction.")
  if (any(!is.finite(road_nrd_year$theta))) stop("Non-finite theta in NRD wind extraction.")

  # 5.5) Downscale using mean-conserving formula
  nrd1 <- downscale_one_layer_formula(nrd12, pmgrid, idx, road_nrd_year, eps = eps, align_fun = align_fun)
  onr1 <- downscale_one_layer_formula(onr12, pmgrid, idx, road_onr_year, eps = eps, align_fun = align_fun)

  # 5.6) Save NRD / ONR
  f_nrd <- file.path(out_data, sprintf("%s_1km_NRD_%s_%d.fst", state_abbr, metric, yr))
  f_onr <- file.path(out_data, sprintf("%s_1km_ONR_%s_%d.fst", state_abbr, metric, yr))
  save_fst(nrd1, f_nrd)
  save_fst(onr1, f_onr)

  # 5.7) TOTAL = NRD + ONR
  setkey(nrd1, x, y)
  setkey(onr1, x, y)
  tot1 <- nrd1[onr1]
  tot1[, pm1_total := pm1 + i.pm1]
  tot1 <- tot1[, .(x, y, pm1_total)]

  f_tot <- file.path(out_data, sprintf("%s_1km_TOTAL_%s_%d.fst", state_abbr, metric, yr))
  save_fst(tot1, f_tot)

  msg("Ranges (NRD pm1):   ", paste(range(nrd1$pm1, na.rm = TRUE), collapse = " ~ "))
  msg("Ranges (ONR pm1):   ", paste(range(onr1$pm1, na.rm = TRUE), collapse = " ~ "))
  msg("Ranges (TOTAL pm1): ", paste(range(tot1$pm1_total, na.rm = TRUE), collapse = " ~ "))

  rm(nrd12, onr12, mm, pmgrid, idx, wind_nrd, wind_onr,
     road_nrd_year, road_onr_year, nrd1, onr1, tot1)
  gc()
}

# 6) Optional: panel maps (only if ALL needed files exist)
if (isTRUE(make_panel_maps)) {
  msg("\n=== MAKING PANEL MAPS (2011–2020) ===")

  needed_nrd <- file.path(out_data, sprintf("%s_1km_NRD_%s_%d.fst", state_abbr, metric, years))
  needed_onr <- file.path(out_data, sprintf("%s_1km_ONR_%s_%d.fst", state_abbr, metric, years))
  needed_tot <- file.path(out_data, sprintf("%s_1km_TOTAL_%s_%d.fst", state_abbr, metric, years))

  if (!all(file.exists(c(needed_nrd, needed_onr, needed_tot)))) {
    msg("Panel maps skipped: some yearly FST files are missing.")
  } else {
    stack_nrd <- read_stack_years("NRD", years, metric, out_data)
    stack_onr <- read_stack_years("ONR", years, metric, out_data)
    stack_tot <- read_stack_years("TOTAL", years, metric, out_data)

    pdf_nrd <- file.path(out_fig, sprintf("PANEL_%s_NRD_%s_2011_2020.pdf", state_abbr, metric))
    pdf_onr <- file.path(out_fig, sprintf("PANEL_%s_ONR_%s_2011_2020.pdf", state_abbr, metric))
    pdf_tot <- file.path(out_fig, sprintf("PANEL_%s_TOTAL_%s_2011_2020.pdf", state_abbr, metric))

    plot_panel_state(stack_nrd, pdf_nrd, state_abbr, p4s, facet_ncol, legend_upper_q)
    plot_panel_state(stack_onr, pdf_onr, state_abbr, p4s, facet_ncol, legend_upper_q)
    plot_panel_state(stack_tot, pdf_tot, state_abbr, p4s, facet_ncol, legend_upper_q)

    rm(stack_nrd, stack_onr, stack_tot); gc()
  }
}

msg("\nALL DONE ✅")
msg("Final data (FST): ", out_data)
msg("Panel maps (PDF): ", out_fig)
