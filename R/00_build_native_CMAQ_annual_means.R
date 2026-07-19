# ==============================================================================
# Build annual and 2011–2020 mean native CMAQ-ISAM PM2.5 surfaces
#
# Raw inputs:
#   PM25_TOT_NRD_cmaq_YYYY-01_YYYY-12.rds
#   PM25_TOT_ONR_cmaq_YYYY-01_YYYY-12.rds
#
# Raw columns:
#   x, y, PM25_TOT_NRD or PM25_TOT_ONR, Date
#
# Outputs:
#   ANNUAL_NRD_mean_YYYY.rds
#   ANNUAL_ONR_mean_YYYY.rds
#   DECADE_NRD_mean_2011_2020.rds
#   DECADE_ONR_mean_2011_2020.rds
#   CMAQ_annual_processing_QC.csv
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------------------------

raw_dir <- paste0(
  "/projects/HAQ_LAB/Sumaiya/cmaq/cmaq_output/",
  "R/data/annual_rds/median_ratio"
)

output_dir <- paste0(
  "/scratch/xshan2/R_Code/Automobiles/",
  "FIGURES/NATIVE_CMAQ_REBUILT"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

years   <- 2011:2020
sectors <- c("NRD", "ONR")

overwrite_existing <- FALSE

# ------------------------------------------------------------------------------
# 2. Helper: locate exactly one raw file
# ------------------------------------------------------------------------------

find_raw_file <- function(sector, year) {

  expected_name <- sprintf(
    "PM25_TOT_%s_cmaq_%d-01_%d-12.rds",
    sector,
    year,
    year
  )

  expected_path <- file.path(
    raw_dir,
    expected_name
  )

  if (!file.exists(expected_path)) {
    stop(
      "Required raw RDS file does not exist:\n",
      expected_path
    )
  }

  expected_path
}

# ------------------------------------------------------------------------------
# 3. Helper: process one sector-year
# ------------------------------------------------------------------------------

process_sector_year <- function(sector, year) {

  raw_file <- find_raw_file(
    sector = sector,
    year = year
  )

  value_column <- paste0(
    "PM25_TOT_",
    sector
  )

  output_file <- file.path(
    output_dir,
    sprintf(
      "ANNUAL_%s_mean_%d.rds",
      sector,
      year
    )
  )

  if (
    file.exists(output_file) &&
    !overwrite_existing
  ) {

    message(
      "Using existing annual file:\n",
      output_file
    )

    annual_dt <- as.data.table(
      readRDS(output_file)
    )

    return(
      list(
        annual_data = annual_dt,
        qc = data.table(
          sector = sector,
          year = year,
          raw_file = raw_file,
          annual_file = output_file,
          status = "existing",
          n_daily_rows = NA_integer_,
          n_dates = NA_integer_,
          n_grid_cells = nrow(annual_dt),
          daily_min = NA_real_,
          daily_mean = NA_real_,
          daily_max = NA_real_,
          annual_min = min(annual_dt$value, na.rm = TRUE),
          annual_mean = mean(annual_dt$value, na.rm = TRUE),
          annual_median = median(annual_dt$value, na.rm = TRUE),
          annual_p95 = as.numeric(
            quantile(
              annual_dt$value,
              0.95,
              na.rm = TRUE
            )
          ),
          annual_p99 = as.numeric(
            quantile(
              annual_dt$value,
              0.99,
              na.rm = TRUE
            )
          ),
          annual_max = max(annual_dt$value, na.rm = TRUE)
        )
      )
    )
  }

  cat(
    "\n========================================\n",
    "PROCESSING ",
    sector,
    " ",
    year,
    "\n",
    "========================================\n",
    "Input:\n",
    raw_file,
    "\n",
    sep = ""
  )

  dt <- readRDS(raw_file)
  setDT(dt)

  required_columns <- c(
    "x",
    "y",
    value_column,
    "Date"
  )

  missing_columns <- setdiff(
    required_columns,
    names(dt)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing required columns in:\n",
      raw_file,
      "\nMissing columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }

  # Keep only required columns to reduce memory
  dt <- dt[, .(
    x = as.numeric(x),
    y = as.numeric(y),
    Date = as.IDate(Date),
    value = as.numeric(
      get(value_column)
    )
  )]

  n_daily_rows <- nrow(dt)
  n_dates <- uniqueN(dt$Date)

  if (n_dates < 360 || n_dates > 366) {
    warning(
      sector,
      " ",
      year,
      " contains ",
      n_dates,
      " unique dates."
    )
  }

  if (
    any(!is.finite(dt$x)) ||
    any(!is.finite(dt$y))
  ) {
    stop(
      "Non-finite x or y coordinates found in:\n",
      raw_file
    )
  }

  # Convert non-finite concentration values to NA
  dt[
    !is.finite(value),
    value := NA_real_
  ]

  # Negative concentrations should not be retained
  n_negative <- dt[
    value < 0,
    .N
  ]

  if (n_negative > 0) {
    warning(
      n_negative,
      " negative values were set to zero in ",
      sector,
      " ",
      year,
      "."
    )

    dt[
      value < 0,
      value := 0
    ]
  }

  daily_min  <- min(dt$value, na.rm = TRUE)
  daily_mean <- mean(dt$value, na.rm = TRUE)
  daily_max  <- max(dt$value, na.rm = TRUE)

  # Annual mean of daily concentrations for each CMAQ grid cell
  annual_dt <- dt[, .(
    value = mean(
      value,
      na.rm = TRUE
    ),
    n_days = sum(
      is.finite(value)
    )
  ), by = .(
    x,
    y
  )]

  annual_dt[
    !is.finite(value),
    value := NA_real_
  ]

  if (anyDuplicated(
    annual_dt[, .(
      x,
      y
    )]
  ) > 0) {
    stop(
      "Duplicate grid coordinates remain after annual aggregation for ",
      sector,
      " ",
      year,
      "."
    )
  }

  setorder(
    annual_dt,
    x,
    y
  )

  saveRDS(
    annual_dt,
    output_file,
    compress = TRUE
  )

  if (!file.exists(output_file)) {
    stop(
      "Annual output was not created:\n",
      output_file
    )
  }

  qc <- data.table(
    sector = sector,
    year = year,
    raw_file = raw_file,
    annual_file = output_file,
    status = "generated",
    n_daily_rows = n_daily_rows,
    n_dates = n_dates,
    n_grid_cells = nrow(annual_dt),
    minimum_days_per_cell = min(
      annual_dt$n_days,
      na.rm = TRUE
    ),
    maximum_days_per_cell = max(
      annual_dt$n_days,
      na.rm = TRUE
    ),
    daily_min = daily_min,
    daily_mean = daily_mean,
    daily_max = daily_max,
    annual_min = min(
      annual_dt$value,
      na.rm = TRUE
    ),
    annual_mean = mean(
      annual_dt$value,
      na.rm = TRUE
    ),
    annual_median = median(
      annual_dt$value,
      na.rm = TRUE
    ),
    annual_p95 = as.numeric(
      quantile(
        annual_dt$value,
        0.95,
        na.rm = TRUE
      )
    ),
    annual_p99 = as.numeric(
      quantile(
        annual_dt$value,
        0.99,
        na.rm = TRUE
      )
    ),
    annual_max = max(
      annual_dt$value,
      na.rm = TRUE
    )
  )

  cat(
    "Created:\n",
    output_file,
    "\nGrid cells: ",
    format(
      nrow(annual_dt),
      big.mark = ","
    ),
    "\nDates: ",
    n_dates,
    "\n",
    sep = ""
  )

  rm(dt)
  gc()

  list(
    annual_data = annual_dt,
    qc = qc
  )
}

# ------------------------------------------------------------------------------
# 4. Process all files and accumulate annual summaries
# ------------------------------------------------------------------------------

all_qc <- list()
annual_by_sector <- list()

for (sector in sectors) {

  sector_annual_list <- vector(
    mode = "list",
    length = length(years)
  )

  for (i in seq_along(years)) {

    year <- years[i]

    result <- process_sector_year(
      sector = sector,
      year = year
    )

    annual_dt <- copy(
      result$annual_data
    )

    annual_dt[, year := year]
    annual_dt[, sector := sector]

    sector_annual_list[[i]] <- annual_dt

    all_qc[[length(all_qc) + 1L]] <- result$qc

    rm(
      result,
      annual_dt
    )

    gc()
  }

  annual_by_sector[[sector]] <- rbindlist(
    sector_annual_list,
    use.names = TRUE
  )

  rm(sector_annual_list)
  gc()
}

# ------------------------------------------------------------------------------
# 5. Calculate 2011–2020 mean for each sector
# ------------------------------------------------------------------------------

decade_results <- list()

for (sector in sectors) {

  sector_stack <- annual_by_sector[[sector]]

  decade_dt <- sector_stack[, .(
    value = mean(
      value,
      na.rm = TRUE
    ),
    n_years = sum(
      is.finite(value)
    )
  ), by = .(
    x,
    y
  )]

  decade_dt[, sector := sector]

  decade_file <- file.path(
    output_dir,
    sprintf(
      "DECADE_%s_mean_2011_2020.rds",
      sector
    )
  )

  saveRDS(
    decade_dt,
    decade_file,
    compress = TRUE
  )

  decade_results[[sector]] <- decade_dt

  cat(
    "\nCreated decade mean:\n",
    decade_file,
    "\nGrid cells: ",
    format(
      nrow(decade_dt),
      big.mark = ","
    ),
    "\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# 6. Create total on-road PM2.5
# ------------------------------------------------------------------------------

nrd_decade <- copy(
  decade_results[["NRD"]]
)[, .(
  x,
  y,
  NRD = value,
  n_years_NRD = n_years
)]

onr_decade <- copy(
  decade_results[["ONR"]]
)[, .(
  x,
  y,
  ONR = value,
  n_years_ONR = n_years
)]

total_decade <- merge(
  nrd_decade,
  onr_decade,
  by = c(
    "x",
    "y"
  ),
  all = TRUE,
  sort = TRUE
)

if (
  anyNA(total_decade$NRD) ||
  anyNA(total_decade$ONR)
) {
  warning(
    "NRD and ONR do not have identical grid coverage."
  )
}

total_decade[, value :=
               NRD + ONR]

total_decade[, n_years :=
               pmin(
                 n_years_NRD,
                 n_years_ONR
               )]

total_decade <- total_decade[, .(
  x,
  y,
  value,
  n_years,
  sector = "TOTAL"
)]

total_file <- file.path(
  output_dir,
  "DECADE_TOTAL_mean_2011_2020.rds"
)

saveRDS(
  total_decade,
  total_file,
  compress = TRUE
)

decade_results[["TOTAL"]] <- total_decade

# ------------------------------------------------------------------------------
# 7. Save QC
# ------------------------------------------------------------------------------

qc_table <- rbindlist(
  all_qc,
  use.names = TRUE,
  fill = TRUE
)

qc_file <- file.path(
  output_dir,
  "CMAQ_annual_processing_QC.csv"
)

fwrite(
  qc_table,
  qc_file
)

decade_qc <- rbindlist(
  lapply(
    names(decade_results),
    function(sector) {

      dt <- decade_results[[sector]]

      data.table(
        sector = sector,
        n_grid_cells = nrow(dt),
        mean = mean(
          dt$value,
          na.rm = TRUE
        ),
        median = median(
          dt$value,
          na.rm = TRUE
        ),
        p05 = as.numeric(
          quantile(
            dt$value,
            0.05,
            na.rm = TRUE
          )
        ),
        p95 = as.numeric(
          quantile(
            dt$value,
            0.95,
            na.rm = TRUE
          )
        ),
        p99 = as.numeric(
          quantile(
            dt$value,
            0.99,
            na.rm = TRUE
          )
        ),
        maximum = max(
          dt$value,
          na.rm = TRUE
        ),
        minimum_years = min(
          dt$n_years,
          na.rm = TRUE
        ),
        maximum_years = max(
          dt$n_years,
          na.rm = TRUE
        )
      )
    }
  ),
  use.names = TRUE
)

decade_qc_file <- file.path(
  output_dir,
  "CMAQ_decade_mean_QC.csv"
)

fwrite(
  decade_qc,
  decade_qc_file
)

print(qc_table)
print(decade_qc)

cat(
  "\n========================================\n",
  "PROCESSING COMPLETED SUCCESSFULLY\n",
  "========================================\n",
  "Output directory:\n",
  output_dir,
  "\n\nAnnual QC:\n",
  qc_file,
  "\n\nDecade QC:\n",
  decade_qc_file,
  "\n",
  sep = ""
)
