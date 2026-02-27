#!/usr/bin/env Rscript
# ============================================================
# Full pipeline (ALL IN ONE):
#  (1) Mean-conserving 1-km downscaling of 12-km CMAQ annual PM2.5
#      using roadiness + wind + wind-road alignment + AADT^alpha
#  (2) Panel plotting (2011–2020) for NRD / ONR / TOTAL
#
# Key formula (for each 1-km cell i within CMAQ grid cell g):
#
#   Atilde_i = (A_i + eps) / (mean_{j in g}(A_j) + eps)   # dimensionless within-grid AADT
#   S_i      = r_i * U_i^{-1} * f(Δθ_i) * (Atilde_i^alpha)
#   C_i      = C_g * S_i / mean_{j in g}(S_j)             # mean-conserving redistribution
#
# Improvements included for your plotting request:
#   - Less “all purple” look: use a tighter upper quantile (e.g., 0.98) for the color limit
#   - Smaller legend: control colorbar height/width and fewer breaks
#
# Output:
#   Downscaled FSTs:
#     DATA_FST/VA_1km_{NRD,ONR,TOTAL}_{metric}_{year}_a{alpha}.fst
#   Panel PDFs:
#     FIGURES/PANEL_{state}_{layer}_{metric}_2011_2020_a{alpha}.pdf
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
  library(grid)      # for unit()
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
annual_dir <- file.path(base_dir, "FIGURES")  # ANNUAL_{NRD,ONR}_{mean|median}_{year}.rds

# ---- 1 km roadiness grid (must contain x,y and roadiness columns)
roadiness_fst <- sprintf("/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/%s/roadiness_1km_hw_loc_%s.fst",
                         state_abbr, state_abbr)

# ---- Roads for computing phi (bearing)
roads_shp <- "/scratch/xshan2/R_Code/Automobiles/roadiness/TRAN_Virginia_State_Shape/Shape/Trans_RoadSegment.shp"
roads_filter_expr <- NULL

# ---- AADT shapefile folder (unzipped)
aadt_base_dir <- "/scratch/xshan2/R_Code/Automobiles/AADT"
aadt_years_available <- 2011:2017

# For missing AADT years (2018–2020):
# "carry_last" -> reuse 2017 AADT
# "ones"       -> set AADT=1 so Atilde=1 (no AADT effect)
aadt_missing_strategy <- "carry_last"

# ---- Wind GeoTIFF: annual mean u/v (Band1=u, Band2=v)
wind_dir     <- sprintf("/scratch/xshan2/NLDAS_wind/%s", state_abbr)
wind_tif_fmt <- file.path(wind_dir, sprintf("NLDAS_%s_%%d_Wind_EN_mean.tif", state_abbr))

# ---- Years and CMAQ metric
years  <- 2011:2020
metric <- "mean"  # or "median"

# ---- Roadiness columns
road_col_NRD <- "nroad.dist2_loc"
road_col_ONR <- "nroad.dist2_hw"

# ---- Alignment function choice
align_fun <- "cos01"  # "cos01" or "cospos"

# ---- Chunk size
chunk_size <- 50000

# ---- Numerical safeguard
eps <- 1e-12

# ---- Alpha sensitivity values
alpha_values <- c(0.5, 1, 2)

# ---- Output folders
out_dir  <- file.path(base_dir, sprintf("DOWNSCALE_1KM_%s_ROAD_WINDDIR_AADT_NORM_ALPHA", state_abbr))
out_data <- file.path(out_dir, "DATA_FST")
out_fig  <- file.path(out_dir, "FIGURES")
out_aux  <- file.path(out_dir, "AUX")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_aux,  recursive = TRUE, showWarnings = FALSE)

phi_cache_fst  <- file.path(out_aux, sprintf("%s_1km_phi_from_roads.fst", state_abbr))
aadt_cache_dir <- file.path(out_aux, "AADT_CACHE")
dir.create(aadt_cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Plot controls (YOUR REQUEST: smaller legend + less purple)
make_panel_maps <- TRUE
panel_alphas_to_plot <- c(1)   # usually alpha=1 for main result; set c(0.5,1,2) if you want all
facet_ncol <- 5

# IMPORTANT: lowering this makes the map less “all purple”
# 0.98 is a good start; 0.95 is more contrast; 0.99 is closer to what you had.
legend_upper_q <- 0.98

# Make legend smaller
legend_bar_height_mm <- 28   # smaller than default
legend_bar_width_mm  <- 5
legend_n_breaks <- 4

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
# Clip 1 km points to state boundary
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
# Load 1 km roadiness grid (x,y,road)
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
# Build unique CMAQ centroid grid and map 1-km points -> nearest CMAQ centroid
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
# Bearing (phi) from LINESTRING geometry (radians)
# ------------------------------------------------------------
bearing_from_linestring <- function(geom) {
  coords <- st_coordinates(geom)
  if (nrow(coords) < 2) return(NA_real_)
  x1 <- coords[1, "X"];  y1 <- coords[1, "Y"]
  x2 <- coords[nrow(coords), "X"]; y2 <- coords[nrow(coords), "Y"]
  atan2(y2 - y1, x2 - x1)
}

# ------------------------------------------------------------
# Compute or load phi (bearing of nearest road segment) for each 1-km point
# ------------------------------------------------------------
compute_or_load_phi <- function(road_xy_dt, roads_shp, p4s, phi_cache_fst,
                                chunk_size = 50000, filter_expr = NULL) {
  stopif_has_cols(road_xy_dt, c("x","y"), where = "road_xy_dt for phi")

  if (file.exists(phi_cache_fst)) {
    msg("Loading cached phi: ", phi_cache_fst)
    phi_dt <- read_fst(phi_cache_fst, as.data.table = TRUE)
    setDT(phi_dt)
    stopif_has_cols(phi_dt, c("x","y","phi"), where = "phi_cache_fst")
    return(phi_dt)
  }

  if (!file.exists(roads_shp)) stop("roads_shp not found: ", roads_shp)

  msg("Reading roads for phi: ", roads_shp)
  roads <- st_read(roads_shp, quiet = TRUE)

  # Drop Z/M dimensions if present (safe even if none)
  roads <- st_zm(roads, drop = TRUE, what = "ZM")

  if (nrow(roads) == 0) stop("roads_shp has 0 features: ", roads_shp)

  if (!is.null(filter_expr)) {
    roads <- roads[eval(filter_expr, envir = roads), ]
    if (nrow(roads) == 0) stop("No roads left after filtering. Check roads_filter_expr.")
  }

  roads <- st_transform(roads, st_crs(p4s))
  roads <- st_cast(roads, "LINESTRING", warn = FALSE)

  msg("Computing road bearings (phi) for ", nrow(roads), " segments ...")
  roads$phi <- vapply(st_geometry(roads), bearing_from_linestring, numeric(1))
  if (sum(!is.finite(roads$phi)) > 0) stop("Some road bearings are NA/Inf. Check road geometry.")

  pts <- st_as_sf(road_xy_dt, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

  n <- nrow(pts)
  idx <- integer(n)

  msg("Computing nearest-road index for phi: n=", n, " (chunk_size=", chunk_size, ") ...")
  for (s in seq(1, n, by = chunk_size)) {
    e <- min(s + chunk_size - 1, n)
    idx[s:e] <- st_nearest_feature(pts[s:e, ], roads)
    msg("  done: ", e, " / ", n)
  }

  phi_dt <- copy(road_xy_dt)
  phi_dt[, phi := as.numeric(roads$phi[idx])]
  if (sum(!is.finite(phi_dt$phi)) > 0) stop("phi_dt contains NA/Inf, unexpected.")

  write_fst(phi_dt[, .(x,y,phi)], phi_cache_fst)
  msg("WROTE phi cache: ", phi_cache_fst, " | rows=", nrow(phi_dt))

  phi_dt[, .(x,y,phi)]
}

# ------------------------------------------------------------
# Extract annual mean wind u/v -> U and theta at 1-km points
# (Transform points to wind raster CRS)
# ------------------------------------------------------------
extract_wind_uv_to_road <- function(road_xy_dt, yr, p4s, eps = 1e-12, method = "bilinear") {
  wind_tif <- sprintf(wind_tif_fmt, yr)
  if (!file.exists(wind_tif)) stop("Missing wind tif: ", wind_tif)

  wd <- raster::brick(wind_tif)
  if (raster::nlayers(wd) < 2) stop("Wind tif must have >=2 bands (u,v): ", wind_tif)

  rcrs_txt <- as.character(raster::crs(wd))
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
    stop("Wind extraction has non-finite u/v for ", sum(bad), " points. ",
         "Check wind tif mask/coverage and state clipping.")
  }

  out[, `:=`(
    U     = sqrt(u^2 + v^2) + eps,
    theta = atan2(v, u)
  )]
  out
}

# ------------------------------------------------------------
# Alignment factor f(Δθ)
# ------------------------------------------------------------
align_factor <- function(delta, mode = c("cos01","cospos")) {
  mode <- match.arg(mode)
  if (mode == "cos01")  return((1 + cos(delta)) / 2)
  if (mode == "cospos") return(pmax(cos(delta), 0))
}

# ------------------------------------------------------------
# AADT helpers
# ------------------------------------------------------------
find_shp_in_year_dir <- function(year_dir) {
  shp_files <- list.files(year_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
  if (length(shp_files) == 0) return(NA_character_)
  yr <- sub(".*?(\\d{4}).*", "\\1", basename(year_dir))
  cand <- shp_files[grepl(yr, basename(shp_files), ignore.case = TRUE)]
  if (length(cand) > 0) return(cand[1])
  shp_files[1]
}

pick_aadt_field <- function(nm) {
  nm_low <- tolower(nm)

  # exact "aadt" first
  idx_exact <- which(nm_low == "aadt")
  if (length(idx_exact) > 0) return(nm[idx_exact[1]])

  # then any field containing "aadt", excluding "truck"
  idx_contains <- which(grepl("aadt", nm_low) & !grepl("truck", nm_low))
  if (length(idx_contains) > 0) return(nm[idx_contains[1]])

  NA_character_
}

read_aadt_roads_year <- function(yr, aadt_base_dir, p4s) {
  year_dir <- file.path(aadt_base_dir, paste0("virginia", yr))
  shp <- find_shp_in_year_dir(year_dir)
  if (is.na(shp) || !file.exists(shp)) return(NULL)

  msg("Reading AADT roads: ", shp)
  rd <- st_read(shp, quiet = TRUE)

  # Drop Z/M dimensions (GEOS cannot handle XYM/XYZM)
  rd <- st_zm(rd, drop = TRUE, what = "ZM")
  rd <- rd[!st_is_empty(rd), ]

  rd <- st_transform(rd, st_crs(p4s))
  rd <- st_cast(rd, "LINESTRING", warn = FALSE)

  aadt_field <- pick_aadt_field(names(rd))
  if (is.na(aadt_field)) {
    msg("Available fields:")
    msg(paste(names(rd), collapse = ", "))
    stop("No AADT-like field found in AADT shp: ", shp)
  }

  msg("Using AADT field: ", aadt_field)
  rd$AADT <- suppressWarnings(as.numeric(rd[[aadt_field]]))
  rd$AADT[!is.finite(rd$AADT) | is.na(rd$AADT)] <- 0
  rd
}

attach_aadt_nearest <- function(road_xy_dt, yr, aadt_base_dir, p4s,
                                chunk_size = 50000, cache_dir = NULL) {
  stopif_has_cols(road_xy_dt, c("x","y"), where = "road_xy_dt for AADT")

  cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, sprintf("AADT_nearest_%d.fst", yr)) else NULL
  if (!is.null(cache_file) && file.exists(cache_file)) {
    msg("Loading cached nearest AADT: ", cache_file)
    dtc <- read_fst(cache_file, as.data.table = TRUE)
    setDT(dtc)
    stopif_has_cols(dtc, c("x","y","AADT"), where = "AADT cache")
    return(dtc)
  }

  roads_aadt <- read_aadt_roads_year(yr, aadt_base_dir, p4s)
  if (is.null(roads_aadt) || nrow(roads_aadt) == 0) return(NULL)

  pts <- st_as_sf(road_xy_dt, coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

  n <- nrow(pts)
  idx <- integer(n)

  msg("Computing nearest AADT-road index: n=", n, " (chunk_size=", chunk_size, ") ...")
  for (s in seq(1, n, by = chunk_size)) {
    e <- min(s + chunk_size - 1, n)
    idx[s:e] <- st_nearest_feature(pts[s:e, ], roads_aadt)
    msg("  done: ", e, " / ", n)
  }

  out <- copy(road_xy_dt)
  out[, AADT := as.numeric(roads_aadt$AADT[idx])]
  out[!is.finite(AADT) | is.na(AADT), AADT := 0]

  if (!is.null(cache_file)) {
    write_fst(out[, .(x,y,AADT)], cache_file)
    msg("WROTE AADT cache: ", cache_file, " | rows=", nrow(out))
  }

  out[, .(x,y,AADT)]
}

# ------------------------------------------------------------
# Downscale one layer (NRD or ONR) with AADT^alpha
# ------------------------------------------------------------
downscale_one_layer <- function(pm_dt, pmgrid, idx_cmaq, road_dt,
                                eps = 1e-12, align_fun = "cos01", alpha = 1) {
  stopif_has_cols(pm_dt, c("lon","lat","val"), where = "pm_dt")
  stopif_has_cols(pmgrid, c("lon","lat","grid_id"), where = "pmgrid")
  stopif_has_cols(road_dt, c("x","y","road","phi","U","theta","AADT"), where = "road_dt")

  pm_dt2 <- merge(pm_dt[, .(lon, lat, Cg = val)], pmgrid, by = c("lon","lat"), all.x = TRUE)
  if (any(is.na(pm_dt2$grid_id))) stop("Some CMAQ cells could not be matched to pmgrid by lon/lat.")
  setkey(pm_dt2, grid_id)

  dt <- copy(road_dt)
  dt[, grid_id := pmgrid$grid_id[idx_cmaq]]
  setkey(dt, grid_id)

  dt <- pm_dt2[dt]

  dt[, r := pmax(road, 0)]
  dt[, delta := atan2(sin(theta - phi), cos(theta - phi))]
  dt[, f := align_factor(delta, mode = align_fun)]

  # Dimensionless AADT within CMAQ grid cell
  dt[, Abar := mean(AADT, na.rm = TRUE), by = grid_id]
  dt[, Atilde := (AADT + eps) / (Abar + eps)]

  dt[, S := r * (1 / U) * f * (Atilde ^ alpha)]
  dt[, Sbar := mean(S, na.rm = TRUE), by = grid_id]
  dt[, Ci := Cg * (S / (Sbar + eps))]

  dt[, .(x, y, pm1 = Ci)]
}

save_fst_dt <- function(dt, outfile) {
  write_fst(dt, outfile)
  msg("WROTE: ", outfile, " | rows=", nrow(dt), " | size=", file.info(outfile)$size)
}

# ============================================================
# PLOTTING HELPERS
# ============================================================

read_one_fst <- function(path, layer) {
  if (!file.exists(path)) stop("Missing file: ", path)
  dt <- read_fst(path, as.data.table = TRUE)
  setDT(dt)

  if (layer == "TOTAL") {
    if (!("pm1_total" %in% names(dt))) stop("TOTAL file missing pm1_total: ", path)
    setnames(dt, "pm1_total", "val")
  } else {
    if (!("pm1" %in% names(dt))) stop(layer, " file missing pm1: ", path)
    setnames(dt, "pm1", "val")
  }
  dt[, .(x, y, val)]
}

fst_path <- function(layer = c("NRD","ONR","TOTAL"), year, alpha, metric, out_data, state_abbr) {
  layer <- match.arg(layer)
  tag <- sprintf("a%g", alpha)
  fname <- sprintf("%s_1km_%s_%s_%d_%s.fst", state_abbr, layer, metric, year, tag)
  file.path(out_data, fname)
}

read_stack_years_alpha <- function(layer, years, alpha, metric, out_data, state_abbr) {
  lst <- lapply(years, function(yr) {
    f <- fst_path(layer, yr, alpha, metric, out_data, state_abbr)
    dt <- read_one_fst(f, layer)
    dt[, year := yr]
    dt
  })
  rbindlist(lst, use.names = TRUE)
}

get_state_outline <- function(state_abbr, p4s) {
  st_outline <- USAboundaries::us_states(resolution = "low")
  st_outline <- st_outline[st_outline$state_abbr == state_abbr, ]
  if (nrow(st_outline) == 0) stop("State not found: ", state_abbr)
  st_transform(st_as_sf(st_outline), st_crs(p4s))
}

plot_panel_years <- function(stack_dt, state_outline, out_pdf,
                             facet_ncol = 5,
                             legend_upper_q = 0.98,
                             legend_bar_height_mm = 28,
                             legend_bar_width_mm = 5,
                             legend_n_breaks = 4) {

  stack_dt <- stack_dt[is.finite(val) & !is.na(val)]
  stack_dt[, year_f := factor(year, levels = sort(unique(year)))]

  xlim <- range(stack_dt$x, na.rm = TRUE)
  ylim <- range(stack_dt$y, na.rm = TRUE)

  # Tighter upper limit -> better contrast (less "all purple")
  lim_upper <- as.numeric(quantile(stack_dt$val, legend_upper_q, na.rm = TRUE))
  if (!is.finite(lim_upper) || lim_upper <= 0) lim_upper <- max(stack_dt$val, na.rm = TRUE)

  # Fewer breaks
  brks <- pretty(c(0, lim_upper), n = legend_n_breaks)

  p <- ggplot() +
    geom_raster(data = stack_dt, aes(x = x, y = y, fill = val)) +
    geom_sf(data = state_outline, fill = NA, color = "grey35", linewidth = 0.2) +
    coord_sf(crs = st_crs(state_outline), xlim = xlim, ylim = ylim, expand = FALSE) +
    facet_wrap(~ year_f, ncol = facet_ncol) +
    scale_fill_viridis_c(
      name = expression(PM[2.5]~(mu*g/m^3)),
      limits = c(0, lim_upper),
      breaks = brks,
      oob = scales::squish,
      na.value = "transparent",
      guide = guide_colorbar(
        barheight = unit(legend_bar_height_mm, "mm"),
        barwidth  = unit(legend_bar_width_mm,  "mm")
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text  = element_text(size = 10)
    )

  grDevices::pdf(out_pdf, width = 16, height = 6.5, useDingbats = FALSE)
  print(p)
  grDevices::dev.off()
  msg("WROTE: ", out_pdf)
}

# ============================================================
# RUN: DOWNSCALE
# ============================================================

msg("=== ", state_abbr, " | Downscaling + plotting (ALL IN ONE) ===")
msg("Output root: ", out_dir)

# 1) Load roadiness grids
road_nrd <- load_roadiness(road_col_NRD)
road_onr <- load_roadiness(road_col_ONR)

# 1.1) Clip to state boundary
road_nrd <- clip_points_to_state(road_nrd, state_abbr, p4s)
road_onr <- clip_points_to_state(road_onr, state_abbr, p4s)

# 1.2) Align NRD/ONR to identical (x,y)
setkey(road_nrd, x, y)
setkey(road_onr, x, y)
road_onr <- road_onr[road_nrd]
road_nrd <- road_nrd[road_onr]
if (nrow(road_nrd) == 0) stop("0 points after clipping/alignment. Check CRS and state clip.")

msg("Final 1 km grid points: ", nrow(road_nrd))

# 2) Phi cache
phi_dt <- compute_or_load_phi(
  road_xy_dt    = road_nrd[, .(x,y)],
  roads_shp     = roads_shp,
  p4s           = p4s,
  phi_cache_fst = phi_cache_fst,
  chunk_size    = chunk_size,
  filter_expr   = roads_filter_expr
)

# 3) Attach phi
setkey(phi_dt, x, y)
setkey(road_nrd, x, y)
setkey(road_onr, x, y)
road_nrd <- phi_dt[road_nrd]
road_onr <- phi_dt[road_onr]

# 4) Reference sf points for CMAQ mapping
road_sf_ref <- st_as_sf(road_nrd[, .(x,y)], coords = c("x","y"), crs = st_crs(p4s), remove = FALSE)

# 5) Loop over years and alphas
for (yr in years) {
  msg("\n-----------------------------")
  msg("YEAR: ", yr, " | metric: ", metric)
  msg("-----------------------------")

  # CMAQ 12-km annual
  nrd12 <- read_annual_rds("NRD", metric, yr)
  onr12 <- read_annual_rds("ONR", metric, yr)

  # CMAQ mapping index (based on NRD grid)
  mm     <- make_pmgrid_and_mapidx(nrd12, road_sf_ref, p4s)
  pmgrid <- mm$pmgrid
  idx    <- mm$idx

  # Wind extraction at 1-km points
  wind_nrd <- extract_wind_uv_to_road(road_nrd[, .(x,y)], yr, p4s, eps)
  wind_onr <- extract_wind_uv_to_road(road_onr[, .(x,y)], yr, p4s, eps)

  # Build per-year road tables (roadiness + phi + wind)
  road_nrd_year <- road_nrd[, .(x,y,road,phi)]
  road_onr_year <- road_onr[, .(x,y,road,phi)]
  setkey(road_nrd_year, x, y); setkey(road_onr_year, x, y)
  setkey(wind_nrd, x, y);      setkey(wind_onr, x, y)
  road_nrd_year <- wind_nrd[road_nrd_year]
  road_onr_year <- wind_onr[road_onr_year]

  # Attach AADT with missing-year strategy
  aadt_year_use <- yr
  if (!(yr %in% aadt_years_available)) {
    if (aadt_missing_strategy == "carry_last") {
      aadt_year_use <- max(aadt_years_available)
      msg("WARNING: AADT missing for ", yr, " -> using ", aadt_year_use, " (carry_last).")
    } else if (aadt_missing_strategy == "ones") {
      aadt_year_use <- NA_integer_
      msg("WARNING: AADT missing for ", yr, " -> using AADT=1 (no AADT effect).")
    } else {
      stop("Unknown aadt_missing_strategy: ", aadt_missing_strategy)
    }
  }

  if (is.na(aadt_year_use)) {
    aadt_nrd <- road_nrd[, .(x,y)][, AADT := 1]
    aadt_onr <- road_onr[, .(x,y)][, AADT := 1]
  } else {
    aadt_nrd <- attach_aadt_nearest(road_nrd[, .(x,y)], aadt_year_use, aadt_base_dir, p4s,
                                    chunk_size = chunk_size, cache_dir = aadt_cache_dir)
    aadt_onr <- attach_aadt_nearest(road_onr[, .(x,y)], aadt_year_use, aadt_base_dir, p4s,
                                    chunk_size = chunk_size, cache_dir = aadt_cache_dir)
    if (is.null(aadt_nrd) || is.null(aadt_onr)) {
      stop("Failed to attach AADT for year_use=", aadt_year_use, ". Check AADT folders under: ", aadt_base_dir)
    }
  }

  # Merge AADT into road tables
  setkey(aadt_nrd, x, y); setkey(aadt_onr, x, y)
  setkey(road_nrd_year, x, y); setkey(road_onr_year, x, y)
  road_nrd_year <- aadt_nrd[road_nrd_year]
  road_onr_year <- aadt_onr[road_onr_year]

  # Alpha sensitivity loop
  for (a in alpha_values) {
    msg("  Downscaling with alpha = ", a)

    nrd1 <- downscale_one_layer(nrd12, pmgrid, idx, road_nrd_year, eps = eps, align_fun = align_fun, alpha = a)
    onr1 <- downscale_one_layer(onr12, pmgrid, idx, road_onr_year, eps = eps, align_fun = align_fun, alpha = a)

    tag <- sprintf("a%g", a)
    f_nrd <- file.path(out_data, sprintf("%s_1km_NRD_%s_%d_%s.fst", state_abbr, metric, yr, tag))
    f_onr <- file.path(out_data, sprintf("%s_1km_ONR_%s_%d_%s.fst", state_abbr, metric, yr, tag))
    save_fst_dt(nrd1, f_nrd)
    save_fst_dt(onr1, f_onr)

    setkey(nrd1, x, y); setkey(onr1, x, y)
    tot1 <- nrd1[onr1]
    tot1[, pm1_total := pm1 + i.pm1]
    tot1 <- tot1[, .(x, y, pm1_total)]

    f_tot <- file.path(out_data, sprintf("%s_1km_TOTAL_%s_%d_%s.fst", state_abbr, metric, yr, tag))
    save_fst_dt(tot1, f_tot)

    rm(nrd1, onr1, tot1); gc()
  }

  rm(nrd12, onr12, mm, pmgrid, idx, wind_nrd, wind_onr,
     road_nrd_year, road_onr_year, aadt_nrd, aadt_onr)
  gc()
}

msg("\nDownscaling finished ✅")
msg("Final data (FST): ", out_data)

# ============================================================
# RUN: PANEL PLOTS
# ============================================================

if (isTRUE(make_panel_maps)) {
  msg("\n=== MAKING PANEL MAPS (2011–2020) ===")

  state_outline <- get_state_outline(state_abbr, p4s)

  for (a in panel_alphas_to_plot) {
    for (layer in c("NRD","ONR","TOTAL")) {

      # Check all needed files exist
      need <- vapply(years, function(yr) {
        file.exists(fst_path(layer, yr, a, metric, out_data, state_abbr))
      }, logical(1))

      if (!all(need)) {
        msg("Skipping panel for ", layer, " alpha=", a, " (missing some yearly files).")
        next
      }

      msg("Panel plot: ", layer, " | alpha=", a)

      stack_dt <- read_stack_years_alpha(layer, years, a, metric, out_data, state_abbr)

      out_pdf <- file.path(out_fig, sprintf("PANEL_%s_%s_%s_2011_2020_a%g.pdf",
                                            state_abbr, layer, metric, a))

      plot_panel_years(
        stack_dt, state_outline, out_pdf,
        facet_ncol = facet_ncol,
        legend_upper_q = legend_upper_q,
        legend_bar_height_mm = legend_bar_height_mm,
        legend_bar_width_mm  = legend_bar_width_mm,
        legend_n_breaks = legend_n_breaks
      )

      rm(stack_dt); gc()
    }
  }

  msg("Panel maps folder: ", out_fig)
}

msg("\nALL DONE ✅")
