#!/usr/bin/env Rscript
# Data source: https://datahub.transportation.gov/stories/s/3uu4-47sa#hpms-state-level-geospatial-data-(2011---2024)
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(sf)
  library(fst)
})

state_name <- "Virginia"
years <- 2011:2020

out_dir <- "/scratch/xshan2/R_Code/Automobiles/HPMS_AADT"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

page_size <- 2000

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

get_base_url <- function(yr) {
  sprintf(
    "https://geo.dot.gov/server/rest/services/Hosted/%s_%d_PR/FeatureServer/0/query",
    state_name, yr
  )
}

get_feature_count <- function(base_url) {
  url <- paste0(base_url, "?where=1%3D1&returnCountOnly=true&f=json")
  x <- jsonlite::fromJSON(url)
  as.integer(x$count)
}

get_page_geojson <- function(base_url, offset, n = 2000) {
  offset_txt <- sprintf("%.0f", offset)

  url <- paste0(
    base_url,
    "?where=1%3D1",
    "&outFields=objectid,year_record,aadt,aadt_combination,aadt_single_unit,f_system,facility_type,county_code,route_id,begin_point,end_point",
    "&returnGeometry=true",
    "&outSR=4326",
    "&f=geojson",
    "&resultOffset=", offset_txt,
    "&resultRecordCount=", n
  )

  txt <- paste(readLines(url, warn = FALSE), collapse = "\n")
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

download_one_year <- function(yr) {
  msg("==========================================")
  msg("Downloading HPMS AADT:", state_name, yr)
  msg("==========================================")

  base_url <- get_base_url(yr)

  out_rds <- file.path(out_dir, sprintf("%s_%d_HPMS_AADT_sf.rds", state_name, yr))
  out_fst <- file.path(out_dir, sprintf("%s_%d_HPMS_AADT_attributes.fst", state_name, yr))

  if (file.exists(out_rds) && file.info(out_rds)$size > 1000) {
    msg("Already exists, skip:", out_rds)
    return(invisible(NULL))
  }

  tmp_dir <- file.path(out_dir, sprintf("tmp_%s_%d_geojson", state_name, yr))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  n_total <- get_feature_count(base_url)
  msg("Total features:", n_total)

  offsets <- seq(0, n_total - 1, by = page_size)

  for (offset in offsets) {
    tmp_file <- file.path(tmp_dir, sprintf("page_%08.0f.geojson", offset))

    if (file.exists(tmp_file) && file.info(tmp_file)$size > 100) {
      msg("Skip existing page:", offset)
      next
    }

    msg("Downloading offset:", offset, "/", n_total)

    ok <- FALSE

    for (try_i in 1:4) {
      x <- try(get_page_geojson(base_url, offset, page_size), silent = TRUE)

      if (!inherits(x, "try-error") && length(x$features) > 0) {
        writeLines(jsonlite::toJSON(x, auto_unbox = TRUE), tmp_file)
        ok <- TRUE
        break
      }

      msg("Retry", try_i, "failed at offset:", offset)
      Sys.sleep(5)
    }

    if (!ok) stop("Failed after retries at offset: ", offset)
  }

  msg("Reading GeoJSON pages ...")

  files <- sort(list.files(tmp_dir, pattern = "\\.geojson$", full.names = TRUE))

  sf_list <- lapply(files, function(f) {
    st_read(f, quiet = TRUE)
  })

  hpms <- do.call(rbind, sf_list)

  hpms$aadt_total <- suppressWarnings(as.numeric(hpms$aadt))
  hpms$aadt_combination <- suppressWarnings(as.numeric(hpms$aadt_combination))
  hpms$aadt_single_unit <- suppressWarnings(as.numeric(hpms$aadt_single_unit))
  hpms$aadt_truck <- hpms$aadt_combination + hpms$aadt_single_unit

  hpms$aadt_total[!is.finite(hpms$aadt_total) | is.na(hpms$aadt_total)] <- 0
  hpms$aadt_truck[!is.finite(hpms$aadt_truck) | is.na(hpms$aadt_truck)] <- 0

  msg("Rows:", nrow(hpms))
  msg("AADT total:")
  print(summary(hpms$aadt_total))
  msg("AADT truck:")
  print(summary(hpms$aadt_truck))

  saveRDS(hpms, out_rds)

  dt_attr <- as.data.table(st_drop_geometry(hpms))
  fst::write_fst(dt_attr, out_fst)

  msg("WROTE:", out_rds)
  msg("WROTE:", out_fst)
}

for (yr in years) {
  download_one_year(yr)
}

msg("ALL DONE")
