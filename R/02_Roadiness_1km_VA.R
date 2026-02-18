#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(sp)
  library(fst)
  library(fasterize)
  library(raster)
  library(areal)
  library(dplyr)
  library(pbmcapply)
})

## =========================================================
## Settings
## =========================================================

p4s <- "+proj=lcc +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97 +a=6370000 +b=6370000"

state_abbr_use <- "VA"
states_use <- c("Virginia_State")

TNMFRC_hw  <- 1:2
TNMFRC_loc <- 4

out_dir  <- "/scratch/xshan2/R_Code/disperseR/Auto/roadiness_2017/VA"
ckpt_dir <- file.path(out_dir, "_ckpt_VA")
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

# optional: safer data.table behavior (less aggressive optimization)
options(datatable.optimize = 1)

## =========================================================
## Checkpoint helper
## =========================================================

cacheRDS <- function(name, expr, force = FALSE) {
  f <- file.path(ckpt_dir, paste0(name, ".rds"))
  if (!force && file.exists(f)) {
    message("[CKPT] load: ", f)
    return(readRDS(f))
  }
  message("[CKPT] run : ", name)
  obj <- eval.parent(substitute(expr))
  saveRDS(obj, f)
  message("[CKPT] save: ", f)
  obj
}

## =========================================================
## Build county boundary (VA) and templates
## =========================================================

counties_use.sf <- cacheRDS("counties_use_sf", {
  counties.us <- USAboundaries::us_counties()
  counties.us$state_name <- NULL

  counties.us %>%
    st_as_sf() %>%
    filter(state_abbr == state_abbr_use) %>%
    st_transform(crs = st_crs(p4s))
})

# Create template rasters (fast; always rebuild to keep them in memory)
road_grid_1km <- raster(extent(counties_use.sf), res = 1000, crs = p4s, vals = NA)

values(road_grid_1km) <- 0


## =========================================================
## Raster grid polygons (slow) -> cache
## =========================================================

road_grid_1km.sf <- cacheRDS("road_grid_1km_sf_poly", {
  rg <- st_as_sf(rasterToPolygons(road_grid_1km))
  rg$id <- 1:nrow(rg)
  rg$road_count <- 0
  rg$road_leng <- 0
  rg$road_leng_m <- 0
  rg
})

## =========================================================
## Road counter function (with robust numeric conversion)
## =========================================================

roadcounter_fn <- function(states_use, TNMFRC, road_grid.sf) {

  road_grid.sf_out <- data.table::copy(road_grid.sf)

  for (state in states_use) {
    message("[roadcounter] state: ", state)

    names.db2 <- list.files("/scratch/xshan2/R_Code/Automobiles/roadiness/",
                            pattern = state, full.names = TRUE)
    names.db2 <- sort(unique(gsub(".zip|_42|_48|_6", "", names.db2)))

    files2 <- list.files(paste0(names.db2, "/Shape/"),
                         pattern = "\\.shp$", full.names = TRUE)

    for (f in files2) {

      roads <- sf::st_read(dsn = f, quiet = TRUE)
      names(roads)[1:min(26, ncol(roads))] <- toupper(names(roads)[1:min(26, ncol(roads))])

      roads.use <- roads[roads$TNMFRC %in% TNMFRC, ]
      if (nrow(roads.use) == 0) next

      roads.use <- sf::st_transform(roads.use, p4s)

      # meters
      roads.use$shape_leng_m <- as.vector(sf::st_length(roads.use))

      # spatial join: assign each road segment to grid id
      ras.roads <- st_join(roads.use, road_grid.sf_out, join = st_intersects) %>% data.table()

      # normalize names
      names(ras.roads) <- tolower(names(ras.roads))

      # robust numeric conversion for shape_leng (may be character/factor)
      if ("shape_leng" %in% names(ras.roads)) {
        ras.roads[, shape_leng := as.numeric(gsub(",", "", as.character(shape_leng)))]
      } else {
        ras.roads[, shape_leng := NA_real_]
      }

      ras.summary <- ras.roads[, .(
        road_count  = .N,
        road_leng   = sum(shape_leng,   na.rm = TRUE),
        road_leng_m = sum(shape_leng_m, na.rm = TRUE)
      ), by = id]

      tmp <- merge(data.table(road_grid.sf_out[, "id"]),
                   ras.summary, by = "id", all.x = TRUE)
      tmp[is.na(tmp)] <- 0

      road_grid.sf_out$road_count  <- road_grid.sf_out$road_count  + tmp$road_count
      road_grid.sf_out$road_leng   <- road_grid.sf_out$road_leng   + tmp$road_leng
      road_grid.sf_out$road_leng_m <- road_grid.sf_out$road_leng_m + tmp$road_leng_m
    }
  }

  road_grid.sf_out
}

## =========================================================
## Run roadcounter (cache each result)
## =========================================================

road_grid_1km_hw.sf <- cacheRDS("VA_1km_hw_sf", {
  roadcounter_fn(states_use, TNMFRC_hw, road_grid_1km.sf)
})

road_grid_1km_loc.sf <- cacheRDS("VA_1km_loc_sf", {
  roadcounter_fn(states_use, TNMFRC_loc, road_grid_1km.sf)
})

## =========================================================
## Fasterize back to rasters (cache)
## =========================================================

road_grid_1km_hw.r <- cacheRDS("VA_1km_hw_count_r", {
  fasterize(road_grid_1km_hw.sf, road_grid_1km, field = "road_count")
})
leng_grid_1km_hw.r <- cacheRDS("VA_1km_hw_leng_r", {
  fasterize(road_grid_1km_hw.sf, road_grid_1km, field = "road_leng")
})
leng_grid_m_1km_hw.r <- cacheRDS("VA_1km_hw_leng_m_r", {
  fasterize(road_grid_1km_hw.sf, road_grid_1km, field = "road_leng_m")
})

road_grid_1km_loc.r <- cacheRDS("VA_1km_loc_count_r", {
  fasterize(road_grid_1km_loc.sf, road_grid_1km, field = "road_count")
})
leng_grid_1km_loc.r <- cacheRDS("VA_1km_loc_leng_r", {
  fasterize(road_grid_1km_loc.sf, road_grid_1km, field = "road_leng")
})
leng_grid_m_1km_loc.r <- cacheRDS("VA_1km_loc_leng_m_r", {
  fasterize(road_grid_1km_loc.sf, road_grid_1km, field = "road_leng_m")
})

## =========================================================
## Build road grids (dt) (cache)
## =========================================================

road_grid_1km.dt <- cacheRDS("VA_1km_grid_dt", {
  data.table(
    coordinates(road_grid_1km_hw.r),
    road_count4_1km_hw   = values(road_grid_1km_hw.r),
    road_leng4_1km_hw    = values(leng_grid_1km_hw.r),
    road_leng4_m_1km_hw  = values(leng_grid_m_1km_hw.r),
    road_count4_1km_loc  = values(road_grid_1km_loc.r),
    road_leng4_1km_loc   = values(leng_grid_1km_loc.r),
    road_leng4_m_1km_loc = values(leng_grid_m_1km_loc.r)
  )
})

road_grid_1km.dt[is.na(road_grid_1km.dt)] <- 0

## =========================================================
## Roadiness worker (dist^-2)
## =========================================================

roadiness.worker <- function(sf.row, dat.sf, points.sf) {

  site.sf <- dat.sf[sf.row, ]

  site_dist <- as.vector(st_distance(site.sf, points.sf)) / 1000
  site_dist[site_dist < 1] <- 1

  data.table(
    nroad.dist2_hw   = sum(points.sf$road_count4_hw   / site_dist^2),
    leng.dist2_hw    = sum(points.sf$road_leng4_hw    / site_dist^2),
    lengm.dist2_hw   = sum(points.sf$road_leng4_m_hw  / site_dist^2),
    nroad.dist2_loc  = sum(points.sf$road_count4_loc  / site_dist^2),
    leng.dist2_loc   = sum(points.sf$road_leng4_loc   / site_dist^2),
    lengm.dist2_loc  = sum(points.sf$road_leng4_m_loc / site_dist^2)
  )
}

## =========================================================
## Convert grids to sf points (avoid rasterFromXYZ round-trip)
## =========================================================

# unify column names (remove _1km suffix)
road_grid_1km.u <- copy(road_grid_1km.dt)
setnames(road_grid_1km.u,
         old = names(road_grid_1km.u),
         new = gsub("_1km", "", names(road_grid_1km.u)))


road_grid_1km.c.sf <- st_as_sf(road_grid_1km.u, coords = c("x", "y"), crs = p4s, remove = FALSE)


# source points: keep only cells with any road metric > 0
road_grid.sf.pos_1km <- road_grid_1km.c.sf[
  road_grid_1km.c.sf$road_count4_hw > 0 |
    road_grid_1km.c.sf$road_leng4_hw > 0 |
    road_grid_1km.c.sf$road_count4_loc > 0 |
    road_grid_1km.c.sf$road_leng4_loc > 0, ]

## =========================================================
## Compute roadiness for all grid cells (cache)
## =========================================================

grid_roadiness_1km.dt <- cacheRDS("VA_roadiness_1km_dt", {
  res <- pbmcapply::pbmclapply(
    1:nrow(road_grid_1km.c.sf),
    roadiness.worker,
    dat.sf = road_grid_1km.c.sf,
    points.sf = road_grid.sf.pos_1km
  )
  rbindlist(res)
})

## =========================================================
## Combine and (optionally) scale/normalize (cache)
## =========================================================

normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

grid_roadiness_full_1km.dt <- cacheRDS("VA_roadiness_full_1km_dt", {
  base_dt <- cbind(road_grid_1km.u[, .(x, y)], grid_roadiness_1km.dt)

  roadiness_cols <- setdiff(names(base_dt), c("x", "y"))
  scaled <- base_dt[, lapply(.SD, scale), .SDcols = roadiness_cols]
  setnames(scaled, names(scaled), paste0(roadiness_cols, "_scale"))

  normed <- base_dt[, lapply(.SD, normalize), .SDcols = roadiness_cols]
  setnames(normed, names(normed), paste0(roadiness_cols, "_norm"))

  cbind(base_dt, scaled, normed)
})


## =========================================================
## Save final outputs
## =========================================================

file.out_1km <- file.path(out_dir, "roadiness_1km_hw_loc_VA.fst")

write_fst(grid_roadiness_full_1km.dt, file.out_1km)


file.out_1km.csv <- file.path(out_dir, "roadiness_1km_hw_loc_VA.csv")

fwrite(grid_roadiness_full_1km.dt, file = file.out_1km.csv)


message("[DONE] All outputs saved in: ", out_dir)
message("[DONE] Checkpoints saved in: ", ckpt_dir)
