#!/usr/bin/env Rscript
# Data source: https://datahub.transportation.gov/stories/s/3uu4-47sa#hpms-state-level-geospatial-data-(2011---2024)
#!/usr/bin/env Rscript

# ==============================================================================
# Download and standardize Virginia HPMS AADT data, 2011-2020 — robust field detection v2
#
# Why this script exists:
#   FHWA/DOT did not use one stable ArcGIS service name for every year.
#
# Sources:
#   2011-2017: archived FHWA ZIP shapefiles
#   2018:      Hosted/Virginia_2018_PR/FeatureServer/0
#   2019:      Hosted/HPMS_Full_VA_2019/FeatureServer/0
#   2020:      Hosted/HPMS_FULL_VA_2020/FeatureServer/0
#
# Outputs for each year:
#   Virginia_YYYY_HPMS_AADT_sf.rds
#     Full line geometry plus standardized fields:
#       aadt_total
#       aadt_combination
#       aadt_single_unit
#       aadt_truck
#
#   Virginia_YYYY_HPMS_AADT_attributes.fst
#     Attribute-only copy for quick QC.
#
# The script:
#   - does not assume the same URL convention for all years;
#   - validates ArcGIS responses before pagination;
#   - automatically selects the line shapefile containing AADT in ZIP archives;
#   - preserves missing truck AADT as NA rather than incorrectly converting it to 0;
#   - writes a source and QC inventory;
#   - stops on a failed year and does not print a false "ALL DONE" message.
#
# Run from a clean Hopper terminal:
#   Rscript Download_Virginia_HPMS_AADT_2011_2020.R \
#     > Download_Virginia_HPMS_AADT_2011_2020.log 2>&1
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(sf)
  library(fst)
})

run_download <- function() {

  # ---------------------------------------------------------------------------
  # 1. Settings
  # ---------------------------------------------------------------------------

  years <- 2011:2020

  out_dir <- "/scratch/xshan2/R_Code/Automobiles/HPMS_AADT"
  raw_dir <- file.path(out_dir, "_RAW_DOWNLOADS")
  tmp_dir <- file.path(out_dir, "_TMP")
  qc_dir  <- file.path(out_dir, "_QC")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir,  recursive = TRUE, showWarnings = FALSE)

  overwrite_existing <- FALSE
  page_size_default <- 2000L
  request_retries <- 5L
  retry_wait_seconds <- 8L

  completion_file <- file.path(
    out_dir,
    "Virginia_HPMS_AADT_2011_2020_DOWNLOAD_COMPLETE.txt"
  )

  source_inventory_file <- file.path(
    qc_dir,
    "Virginia_HPMS_AADT_source_inventory.csv"
  )

  qc_file <- file.path(
    qc_dir,
    "Virginia_HPMS_AADT_QC_2011_2020.csv"
  )

  unlink(completion_file, force = TRUE)

  msg <- function(...) {
    cat(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "|",
      ...,
      "\n"
    )
    flush.console()
  }

  # Explicit source map. Do not replace this with one generic URL formula.
  source_table <- data.table(
    year = years,
    source_type = c(
      rep("zip", 7),
      rep("arcgis", 3)
    ),
    source_url = c(
      sprintf(
        "https://www.fhwa.dot.gov/policyinformation/hpms/shapefiles/virginia%d.zip",
        2011:2017
      ),
      "https://geo.dot.gov/server/rest/services/Hosted/Virginia_2018_PR/FeatureServer/0",
      "https://geo.dot.gov/server/rest/services/Hosted/HPMS_Full_VA_2019/FeatureServer/0",
      "https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_VA_2020/FeatureServer/0"
    )
  )

  fwrite(source_table, source_inventory_file)

  # ---------------------------------------------------------------------------
  # 2. General helpers
  # ---------------------------------------------------------------------------

  safe_numeric <- function(x) {
    suppressWarnings(
      as.numeric(
        gsub(
          ",",
          "",
          trimws(
            as.character(x)
          )
        )
      )
    )
  }

  clean_names <- function(x) {
    x <- tolower(x)
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    make.unique(x, sep = "_")
  }

  read_url_text <- function(url, retries = request_retries) {

    last_error <- NULL

    for (attempt in seq_len(retries)) {

      result <- tryCatch(
        {
          con <- base::url(url, open = "rb")
          raw <- readBin(
            con,
            what = "raw",
            n = 1024L * 1024L * 500L
          )
          close(con)
          rawToChar(raw)
        },
        error = function(e) {
          last_error <<- conditionMessage(e)
          NULL
        }
      )

      if (!is.null(result) && nzchar(result)) {
        return(result)
      }

      msg(
        "URL request attempt ",
        attempt,
        " failed; retrying in ",
        retry_wait_seconds,
        " s. URL: ",
        url
      )

      Sys.sleep(retry_wait_seconds)
    }

    stop(
      "Failed to read URL after ",
      retries,
      " attempts:\n",
      url,
      "\nLast error: ",
      last_error
    )
  }

  download_binary <- function(url, destination) {

    if (
      file.exists(destination) &&
      is.finite(file.info(destination)$size) &&
      file.info(destination)$size > 1000
    ) {
      msg("Using existing download: ", destination)
      return(invisible(destination))
    }

    tmp <- paste0(destination, ".partial")
    unlink(tmp, force = TRUE)

    status <- tryCatch(
      {
        suppressWarnings(
          download.file(
            url = url,
            destfile = tmp,
            mode = "wb",
            quiet = FALSE,
            method = "libcurl"
          )
        )
      },
      error = function(e) {
        msg("libcurl download failed: ", conditionMessage(e))
        NA_integer_
      }
    )

    if (
      is.na(status) ||
      status != 0 ||
      !file.exists(tmp) ||
      file.info(tmp)$size <= 1000
    ) {
      unlink(tmp, force = TRUE)
      stop(
        "Binary download failed:\n",
        url
      )
    }

    if (!file.rename(tmp, destination)) {
      stop(
        "Downloaded file could not be renamed:\n",
        tmp,
        "\nTo:\n",
        destination
      )
    }

    destination
  }

  validate_line_sf <- function(x, context) {

    if (!inherits(x, "sf")) {
      stop(context, " is not an sf object.")
    }

    if (nrow(x) == 0) {
      stop(context, " contains zero features.")
    }

    geometry_types <- unique(
      as.character(
        st_geometry_type(x)
      )
    )

    if (!any(geometry_types %in% c("LINESTRING", "MULTILINESTRING"))) {
      stop(
        context,
        " does not contain line geometry. Types: ",
        paste(geometry_types, collapse = ", ")
      )
    }

    invisible(TRUE)
  }

  # ---------------------------------------------------------------------------
  # 3. Standardize AADT fields
  # ---------------------------------------------------------------------------

  # Score and select a total-AADT field across changing historical schemas.
  #
  # Early HPMS shapefiles occasionally use names such as AADT_2012 or other
  # truncated variants. The previous script accepted only exact "AADT", which
  # caused the valid 2012 line shapefile to be rejected.
  select_total_aadt_field <- function(
      data,
      attribute_names,
      year) {

    candidate_names <- attribute_names[
      grepl(
        "aadt|annual.*average.*daily.*traffic|(^|_)adt($|_)|traffic.*volume",
        attribute_names,
        ignore.case = TRUE,
        perl = TRUE
      )
    ]

    # Do not mistake truck components, percentages, or metadata for total AADT.
    candidate_names <- candidate_names[
      !grepl(
        "comb|single|truck|aadtt|percent|pct|class|factor|route|year_record",
        candidate_names,
        ignore.case = TRUE,
        perl = TRUE
      )
    ]

    if (length(candidate_names) == 0) {
      return(NA_character_)
    }

    field_scores <- vapply(
      candidate_names,
      function(field_name) {

        values <- safe_numeric(
          data[[field_name]]
        )

        n_finite <- sum(
          is.finite(values)
        )

        n_positive <- sum(
          is.finite(values) &
            values > 0
        )

        n_distinct <- data.table::uniqueN(
          values[
            is.finite(values)
          ]
        )

        # A traffic-volume field must contain a meaningful amount of numeric data.
        if (
          n_finite < max(
            10L,
            ceiling(
              0.001 *
                nrow(data)
            )
          ) ||
            n_positive < 5L ||
            n_distinct < 3L
        ) {
          return(-Inf)
        }

        score <- 0

        if (
          identical(
            field_name,
            "aadt"
          )
        ) {
          score <- score + 1000
        }

        if (
          identical(
            field_name,
            "aadt_total"
          )
        ) {
          score <- score + 950
        }

        if (
          grepl(
            paste0(
              "^aadt_?",
              year,
              "$"
            ),
            field_name,
            ignore.case = TRUE
          )
        ) {
          score <- score + 900
        }

        if (
          grepl(
            "^aadt",
            field_name,
            ignore.case = TRUE
          )
        ) {
          score <- score + 700
        }

        if (
          grepl(
            "annual.*average.*daily.*traffic",
            field_name,
            ignore.case = TRUE
          )
        ) {
          score <- score + 650
        }

        if (
          grepl(
            "(^|_)adt($|_)",
            field_name,
            ignore.case = TRUE,
            perl = TRUE
          )
        ) {
          score <- score + 500
        }

        # Prefer fields with broad numeric coverage.
        score <- score +
          100 *
          (
            n_finite /
              max(
                1,
                nrow(data)
              )
          ) +
          log10(
            n_positive +
              1
          ) +
          0.1 *
          log10(
            n_distinct +
              1
          )

        score
      },
      numeric(1)
    )

    if (
      all(
        !is.finite(
          field_scores
        )
      )
    ) {
      return(NA_character_)
    }

    candidate_names[
      which.max(
        field_scores
      )
    ]
  }

  select_component_field <- function(
      data,
      attribute_names,
      component = c(
        "combination",
        "single"
      )) {

    component <- match.arg(
      component
    )

    patterns <- if (
      component ==
        "combination"
    ) {
      c(
        "aadt.*comb",
        "comb.*aadt",
        "^aadt_com",
        "combination.*traffic"
      )
    } else {
      c(
        "aadt.*single",
        "single.*aadt",
        "^aadt_sing",
        "single.*unit.*traffic"
      )
    }

    candidates <- unique(
      unlist(
        lapply(
          patterns,
          function(pattern) {
            attribute_names[
              grepl(
                pattern,
                attribute_names,
                ignore.case = TRUE,
                perl = TRUE
              )
            ]
          }
        ),
        use.names = FALSE
      )
    )

    if (
      length(
        candidates
      ) ==
        0
    ) {
      return(NA_character_)
    }

    scores <- vapply(
      candidates,
      function(field_name) {

        values <- safe_numeric(
          data[[field_name]]
        )

        n_finite <- sum(
          is.finite(values)
        )

        n_positive <- sum(
          is.finite(values) &
            values > 0
        )

        if (
          n_finite <
            10 ||
            n_positive <
            3
        ) {
          return(-Inf)
        }

        100 *
          (
            n_finite /
              max(
                1,
                nrow(data)
              )
          ) +
          log10(
            n_positive +
              1
          )
      },
      numeric(1)
    )

    if (
      all(
        !is.finite(
          scores
        )
      )
    ) {
      return(NA_character_)
    }

    candidates[
      which.max(
        scores
      )
    ]
  }

  standardize_hpms <- function(hpms, year, source_url) {

    validate_line_sf(
      hpms,
      paste0("HPMS ", year)
    )

    hpms <- st_zm(
      hpms,
      drop = TRUE,
      what = "ZM"
    )

    hpms <- hpms[
      !st_is_empty(hpms),
    ]

    old_names <- names(hpms)
    names(hpms) <- clean_names(old_names)

    geometry_column <- attr(hpms, "sf_column")
    attribute_names <- setdiff(
      names(hpms),
      geometry_column
    )

    total_field <- select_total_aadt_field(
      data = hpms,
      attribute_names = attribute_names,
      year = year
    )

    combination_field <- select_component_field(
      data = hpms,
      attribute_names = attribute_names,
      component = "combination"
    )

    single_field <- select_component_field(
      data = hpms,
      attribute_names = attribute_names,
      component = "single"
    )

    msg(
      "Detected fields for ",
      year,
      ": total=",
      total_field,
      "; combination=",
      combination_field,
      "; single_unit=",
      single_field
    )

    if (is.na(total_field)) {
      stop(
        "No total AADT field was identified for ",
        year,
        ". Available fields:\n",
        paste(attribute_names, collapse = ", ")
      )
    }

    hpms$aadt_total <- safe_numeric(
      hpms[[total_field]]
    )

    hpms$aadt_combination <- if (!is.na(combination_field)) {
      safe_numeric(
        hpms[[combination_field]]
      )
    } else {
      NA_real_
    }

    hpms$aadt_single_unit <- if (!is.na(single_field)) {
      safe_numeric(
        hpms[[single_field]]
      )
    } else {
      NA_real_
    }

    # Truck AADT is defined only where both truck component fields are available.
    hpms$aadt_truck <- if (
      !is.na(combination_field) &&
      !is.na(single_field)
    ) {
      hpms$aadt_combination +
        hpms$aadt_single_unit
    } else {
      NA_real_
    }

    # Non-finite total AADT is retained as NA, not converted to zero.
    hpms$aadt_total[
      !is.finite(hpms$aadt_total)
    ] <- NA_real_

    hpms$aadt_combination[
      !is.finite(hpms$aadt_combination)
    ] <- NA_real_

    hpms$aadt_single_unit[
      !is.finite(hpms$aadt_single_unit)
    ] <- NA_real_

    hpms$aadt_truck[
      !is.finite(hpms$aadt_truck)
    ] <- NA_real_

    # Negative traffic values are invalid and are retained as missing.
    hpms$aadt_total[
      hpms$aadt_total < 0
    ] <- NA_real_

    hpms$aadt_combination[
      hpms$aadt_combination < 0
    ] <- NA_real_

    hpms$aadt_single_unit[
      hpms$aadt_single_unit < 0
    ] <- NA_real_

    hpms$aadt_truck[
      hpms$aadt_truck < 0
    ] <- NA_real_

    hpms$hpms_year <- as.integer(year)
    hpms$hpms_source_url <- source_url
    hpms$hpms_total_field <- total_field
    hpms$hpms_combination_field <- ifelse(
      is.na(combination_field),
      NA_character_,
      combination_field
    )
    hpms$hpms_single_unit_field <- ifelse(
      is.na(single_field),
      NA_character_,
      single_field
    )

    hpms
  }

  # ---------------------------------------------------------------------------
  # 4. ZIP archive years: 2011-2017
  # ---------------------------------------------------------------------------

  read_best_aadt_shapefile <- function(
      shapefiles,
      year) {

    if (
      length(
        shapefiles
      ) ==
        0
    ) {
      stop(
        "No shapefiles were found after extracting ",
        year,
        "."
      )
    }

    candidate_objects <- list()
    candidate_summary <- list()

    for (
      shp in shapefiles
    ) {

      msg(
        "Inspecting shapefile: ",
        shp
      )

      obj <- tryCatch(
        st_read(
          shp,
          quiet = TRUE,
          stringsAsFactors = FALSE
        ),
        error = function(e) {
          msg(
            "Could not read ",
            shp,
            ": ",
            conditionMessage(
              e
            )
          )
          NULL
        }
      )

      if (
        is.null(
          obj
        ) ||
          nrow(
            obj
          ) ==
          0
      ) {
        next
      }

      obj <- st_zm(
        obj,
        drop = TRUE,
        what = "ZM"
      )

      geometry_types <- unique(
        as.character(
          st_geometry_type(
            obj
          )
        )
      )

      is_line <- any(
        grepl(
          "LINESTRING",
          geometry_types,
          fixed = TRUE
        )
      )

      # Some legacy files can report a generic GEOMETRY type. Use topological
      # dimension as a fallback.
      if (
        !is_line
      ) {
        dimensions <- suppressWarnings(
          st_dimension(
            st_geometry(
              obj
            )
          )
        )

        is_line <- any(
          dimensions ==
            1,
          na.rm = TRUE
        )
      }

      obj_test <- obj
      names(
        obj_test
      ) <- clean_names(
        names(
          obj_test
        )
      )

      geometry_column <- attr(
        obj_test,
        "sf_column"
      )

      attribute_names <- setdiff(
        names(
          obj_test
        ),
        geometry_column
      )

      total_field <- select_total_aadt_field(
        data = obj_test,
        attribute_names = attribute_names,
        year = year
      )

      msg(
        "  Geometry types: ",
        paste(
          geometry_types,
          collapse = "|"
        )
      )

      msg(
        "  Candidate total-AADT field: ",
        ifelse(
          is.na(
            total_field
          ),
          "<none>",
          total_field
        )
      )

      if (
        is.na(
          total_field
        )
      ) {
        msg(
          "  Available fields: ",
          paste(
            attribute_names,
            collapse = ", "
          )
        )
      }

      if (
        is_line &&
          !is.na(
            total_field
          )
      ) {

        values <- safe_numeric(
          obj_test[[total_field]]
        )

        candidate_objects[[shp]] <- obj

        candidate_summary[[shp]] <- data.table(
          shapefile = shp,
          n_rows = nrow(
            obj
          ),
          n_finite_aadt = sum(
            is.finite(
              values
            )
          ),
          n_positive_aadt = sum(
            is.finite(
              values
            ) &
              values >
              0
          ),
          total_field = total_field
        )
      }

      rm(
        obj_test
      )
    }

    if (
      length(
        candidate_objects
      ) ==
        0
    ) {
      stop(
        "No extracted line shapefile with a usable total-AADT field was found ",
        "for ",
        year,
        ".\nThe updated detector checked names such as AADT, AADT_YYYY, ",
        "and other numeric AADT variants.\nShapefiles checked:\n",
        paste(
          shapefiles,
          collapse = "\n"
        )
      )
    }

    candidate_table <- rbindlist(
      candidate_summary,
      use.names = TRUE,
      fill = TRUE
    )

    setorder(
      candidate_table,
      -n_finite_aadt,
      -n_positive_aadt,
      -n_rows
    )

    selected_name <- candidate_table$shapefile[1]

    msg(
      "Selected shapefile for ",
      year,
      ": ",
      selected_name,
      " (",
      candidate_table$n_rows[1],
      " rows; field ",
      candidate_table$total_field[1],
      "; ",
      candidate_table$n_finite_aadt[1],
      " finite AADT values)"
    )

    candidate_objects[[selected_name]]
  }

  download_zip_year <- function(year, source_url) {

    year_raw_dir <- file.path(
      raw_dir,
      as.character(year)
    )

    extraction_dir <- file.path(
      year_raw_dir,
      "unzipped"
    )

    dir.create(
      year_raw_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    zip_file <- file.path(
      year_raw_dir,
      sprintf(
        "virginia%d.zip",
        year
      )
    )

    download_binary(
      source_url,
      zip_file
    )

    if (
      !dir.exists(extraction_dir) ||
      length(
        list.files(
          extraction_dir,
          recursive = TRUE,
          all.files = FALSE
        )
      ) == 0
    ) {
      unlink(
        extraction_dir,
        recursive = TRUE,
        force = TRUE
      )

      dir.create(
        extraction_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      unzip_status <- tryCatch(
        {
          unzip(
            zip_file,
            exdir = extraction_dir
          )
          TRUE
        },
        error = function(e) {
          msg(
            "Unzip failed for ",
            zip_file,
            ": ",
            conditionMessage(e)
          )
          FALSE
        }
      )

      if (!unzip_status) {
        stop("Could not unzip: ", zip_file)
      }
    }

    shapefiles <- list.files(
      extraction_dir,
      pattern = "\\.shp$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    read_best_aadt_shapefile(
      shapefiles,
      year
    )
  }

  # ---------------------------------------------------------------------------
  # 5. ArcGIS FeatureServer years: 2018-2020
  # ---------------------------------------------------------------------------

  arcgis_json <- function(url) {

    text <- read_url_text(url)

    parsed <- tryCatch(
      fromJSON(
        text,
        simplifyVector = TRUE
      ),
      error = function(e) {
        stop(
          "ArcGIS returned non-JSON content:\n",
          url,
          "\n",
          conditionMessage(e),
          "\nResponse start:\n",
          substr(text, 1, 1000)
        )
      }
    )

    if (!is.null(parsed$error)) {
      stop(
        "ArcGIS error for:\n",
        url,
        "\n",
        paste(
          capture.output(
            str(parsed$error)
          ),
          collapse = "\n"
        )
      )
    }

    parsed
  }

  arcgis_layer_metadata <- function(layer_url) {

    metadata <- arcgis_json(
      paste0(
        layer_url,
        "?f=pjson"
      )
    )

    object_id_field <- metadata$objectIdField

    if (
      is.null(object_id_field) ||
      length(object_id_field) != 1 ||
      is.na(object_id_field) ||
      !nzchar(object_id_field)
    ) {
      object_id_field <- "OBJECTID"
    }

    max_record_count <- suppressWarnings(
      as.integer(
        metadata$maxRecordCount
      )
    )

    if (
      length(max_record_count) != 1 ||
      is.na(max_record_count) ||
      max_record_count < 1
    ) {
      max_record_count <- page_size_default
    }

    list(
      object_id_field = object_id_field,
      max_record_count = min(
        max_record_count,
        page_size_default
      ),
      metadata = metadata
    )
  }

  arcgis_feature_count <- function(layer_url) {

    query_url <- paste0(
      layer_url,
      "/query?",
      "where=1%3D1",
      "&returnCountOnly=true",
      "&f=json"
    )

    result <- arcgis_json(query_url)

    count <- suppressWarnings(
      as.integer(
        result$count
      )
    )

    if (
      length(count) != 1 ||
      is.na(count) ||
      count < 1
    ) {
      stop(
        "ArcGIS feature count was invalid for:\n",
        layer_url,
        "\nParsed count: ",
        paste(count, collapse = ", ")
      )
    }

    count
  }

  build_arcgis_page_url <- function(
      layer_url,
      offset,
      page_size,
      object_id_field) {

    paste0(
      layer_url,
      "/query?",
      "where=1%3D1",
      "&outFields=*",
      "&returnGeometry=true",
      "&outSR=4326",
      "&returnZ=false",
      "&returnM=false",
      "&f=geojson",
      "&orderByFields=",
      URLencode(
        object_id_field,
        reserved = TRUE
      ),
      "&resultOffset=",
      format(
        offset,
        scientific = FALSE,
        trim = TRUE
      ),
      "&resultRecordCount=",
      page_size
    )
  }

  download_arcgis_year <- function(
      year,
      layer_url) {

    metadata <- arcgis_layer_metadata(
      layer_url
    )

    feature_count <- arcgis_feature_count(
      layer_url
    )

    page_size <- metadata$max_record_count

    msg(
      "ArcGIS ",
      year,
      ": ",
      feature_count,
      " features; page size ",
      page_size,
      "; OID field ",
      metadata$object_id_field
    )

    year_tmp_dir <- file.path(
      tmp_dir,
      paste0(
        "arcgis_",
        year
      )
    )

    dir.create(
      year_tmp_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    offsets <- seq(
      0L,
      feature_count - 1L,
      by = page_size
    )

    page_files <- character(
      length(offsets)
    )

    for (i in seq_along(offsets)) {

      offset <- offsets[i]

      page_file <- file.path(
        year_tmp_dir,
        sprintf(
          "page_%08d.geojson",
          offset
        )
      )

      page_files[i] <- page_file

      if (
        file.exists(page_file) &&
        file.info(page_file)$size > 100
      ) {
        msg(
          "Using cached GeoJSON page: ",
          basename(page_file)
        )
        next
      }

      page_url <- build_arcgis_page_url(
        layer_url = layer_url,
        offset = offset,
        page_size = page_size,
        object_id_field = metadata$object_id_field
      )

      msg(
        "Downloading ",
        year,
        " offset ",
        offset,
        " / ",
        feature_count
      )

      page_text <- read_url_text(
        page_url
      )

      page_json <- tryCatch(
        fromJSON(
          page_text,
          simplifyVector = FALSE
        ),
        error = function(e) {
          stop(
            "Could not parse ArcGIS GeoJSON page for ",
            year,
            ", offset ",
            offset,
            ":\n",
            conditionMessage(e),
            "\nResponse start:\n",
            substr(page_text, 1, 1000)
          )
        }
      )

      if (!is.null(page_json$error)) {
        stop(
          "ArcGIS returned an error for ",
          year,
          ", offset ",
          offset,
          ":\n",
          paste(
            capture.output(
              str(page_json$error)
            ),
            collapse = "\n"
          )
        )
      }

      n_features <- length(
        page_json$features
      )

      if (n_features == 0) {
        stop(
          "ArcGIS returned zero features for ",
          year,
          ", offset ",
          offset,
          "."
        )
      }

      writeLines(
        page_text,
        con = page_file,
        useBytes = TRUE
      )
    }

    msg("Reading downloaded GeoJSON pages for ", year)

    page_objects <- lapply(
      page_files,
      function(page_file) {

        object <- st_read(
          page_file,
          quiet = TRUE,
          stringsAsFactors = FALSE
        )

        if (nrow(object) == 0) {
          stop(
            "Downloaded GeoJSON page contains zero features:\n",
            page_file
          )
        }

        object
      }
    )

    hpms <- do.call(
      rbind,
      page_objects
    )

    if (nrow(hpms) != feature_count) {
      warning(
        "Downloaded ",
        nrow(hpms),
        " rows for ",
        year,
        " but ArcGIS count reported ",
        feature_count,
        ". Duplicate OIDs and service pagination will be checked."
      )
    }

    hpms
  }

  # ---------------------------------------------------------------------------
  # 6. Existing-output validation
  # ---------------------------------------------------------------------------

  existing_output_valid <- function(rds_file, year) {

    if (
      !file.exists(rds_file) ||
      !is.finite(file.info(rds_file)$size) ||
      file.info(rds_file)$size <= 1000
    ) {
      return(FALSE)
    }

    object <- tryCatch(
      readRDS(rds_file),
      error = function(e) NULL
    )

    if (is.null(object) || !inherits(object, "sf")) {
      return(FALSE)
    }

    required <- c(
      "aadt_total",
      "aadt_combination",
      "aadt_single_unit",
      "aadt_truck",
      "hpms_year"
    )

    if (!all(required %in% names(object))) {
      return(FALSE)
    }

    if (
      nrow(object) == 0 ||
      !all(
        unique(
          object$hpms_year
        ) ==
          year
      )
    ) {
      return(FALSE)
    }

    TRUE
  }

  # ---------------------------------------------------------------------------
  # 7. Main annual loop
  # ---------------------------------------------------------------------------

  qc_list <- list()

  for (year in years) {

    msg(
      "============================================================"
    )

    msg(
      "Processing Virginia HPMS AADT ",
      year
    )

    msg(
      "============================================================"
    )

    current_year <- year

    source_type <- source_table[
      year == current_year,
      source_type
    ]

    source_url <- source_table[
      year == current_year,
      source_url
    ]

    if (
      length(source_type) != 1 ||
      length(source_url) != 1
    ) {
      stop(
        "Source table did not resolve uniquely for ",
        year,
        "."
      )
    }

    out_rds <- file.path(
      out_dir,
      sprintf(
        "Virginia_%d_HPMS_AADT_sf.rds",
        year
      )
    )

    out_fst <- file.path(
      out_dir,
      sprintf(
        "Virginia_%d_HPMS_AADT_attributes.fst",
        year
      )
    )

    if (
      !overwrite_existing &&
      existing_output_valid(
        out_rds,
        year
      )
    ) {
      msg(
        "Validated existing standardized output; skipping download: ",
        out_rds
      )

      hpms <- readRDS(
        out_rds
      )

    } else {

      hpms_raw <- if (
        source_type ==
          "zip"
      ) {
        download_zip_year(
          year,
          source_url
        )
      } else if (
        source_type ==
          "arcgis"
      ) {
        download_arcgis_year(
          year,
          source_url
        )
      } else {
        stop(
          "Unknown source type for ",
          year,
          ": ",
          source_type
        )
      }

      hpms <- standardize_hpms(
        hpms_raw,
        year = year,
        source_url = source_url
      )

      saveRDS(
        hpms,
        out_rds,
        compress = "xz"
      )

      attributes <- as.data.table(
        st_drop_geometry(
          hpms
        )
      )

      write_fst(
        attributes,
        out_fst
      )

      if (
        !file.exists(out_rds) ||
        file.info(out_rds)$size <= 1000 ||
        !file.exists(out_fst) ||
        file.info(out_fst)$size <= 100
      ) {
        stop(
          "Standardized outputs were not created correctly for ",
          year,
          "."
        )
      }
    }

    finite_total <- hpms$aadt_total[
      is.finite(
        hpms$aadt_total
      )
    ]

    finite_truck <- hpms$aadt_truck[
      is.finite(
        hpms$aadt_truck
      )
    ]

    qc_list[[as.character(year)]] <- data.table(
      year = year,
      source_type = source_type,
      source_url = source_url,
      n_features = nrow(hpms),
      geometry_types = paste(
        sort(
          unique(
            as.character(
              st_geometry_type(hpms)
            )
          )
        ),
        collapse = "|"
      ),
      crs = st_crs(hpms)$input,
      total_field = unique(
        hpms$hpms_total_field
      )[1],
      combination_field = unique(
        hpms$hpms_combination_field
      )[1],
      single_unit_field = unique(
        hpms$hpms_single_unit_field
      )[1],
      n_finite_total_aadt = length(
        finite_total
      ),
      percent_missing_total_aadt =
        100 *
        mean(
          !is.finite(
            hpms$aadt_total
          )
        ),
      total_aadt_median = if (
        length(finite_total) >
          0
      ) {
        median(
          finite_total
        )
      } else {
        NA_real_
      },
      total_aadt_p95 = if (
        length(finite_total) >
          0
      ) {
        as.numeric(
          quantile(
            finite_total,
            0.95
          )
        )
      } else {
        NA_real_
      },
      total_aadt_maximum = if (
        length(finite_total) >
          0
      ) {
        max(
          finite_total
        )
      } else {
        NA_real_
      },
      n_finite_truck_aadt = length(
        finite_truck
      ),
      percent_missing_truck_aadt =
        100 *
        mean(
          !is.finite(
            hpms$aadt_truck
          )
        ),
      truck_aadt_median = if (
        length(finite_truck) >
          0
      ) {
        median(
          finite_truck
        )
      } else {
        NA_real_
      },
      truck_aadt_p95 = if (
        length(finite_truck) >
          0
      ) {
        as.numeric(
          quantile(
            finite_truck,
            0.95
          )
        )
      } else {
        NA_real_
      },
      truck_aadt_maximum = if (
        length(finite_truck) >
          0
      ) {
        max(
          finite_truck
        )
      } else {
        NA_real_
      },
      rds_file = out_rds,
      rds_size_mb = round(
        file.info(out_rds)$size /
          1024^2,
        3
      ),
      fst_file = out_fst,
      fst_size_mb = round(
        file.info(out_fst)$size /
          1024^2,
        3
      )
    )

    msg(
      "Completed ",
      year,
      ": ",
      format(
        nrow(hpms),
        big.mark = ","
      ),
      " line features"
    )

    rm(hpms)
    gc()
  }

  # ---------------------------------------------------------------------------
  # 8. Final completeness checks
  # ---------------------------------------------------------------------------

  qc <- rbindlist(
    qc_list,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    qc,
    year
  )

  fwrite(
    qc,
    qc_file
  )

  if (
    !identical(
      as.integer(
        qc$year
      ),
      as.integer(
        years
      )
    )
  ) {
    stop(
      "Final QC does not contain exactly 2011-2020."
    )
  }

  output_rds <- file.path(
    out_dir,
    sprintf(
      "Virginia_%d_HPMS_AADT_sf.rds",
      years
    )
  )

  output_fst <- file.path(
    out_dir,
    sprintf(
      "Virginia_%d_HPMS_AADT_attributes.fst",
      years
    )
  )

  missing_outputs <- c(
    output_rds,
    output_fst,
    source_inventory_file,
    qc_file
  )

  missing_outputs <- missing_outputs[
    !file.exists(
      missing_outputs
    )
  ]

  if (length(missing_outputs) > 0) {
    stop(
      "Missing required outputs:\n",
      paste(
        missing_outputs,
        collapse = "\n"
      )
    )
  }

  completion_text <- c(
    "Virginia HPMS AADT download and standardization completed.",
    paste0(
      "Years: ",
      paste(
        years,
        collapse = ", "
      )
    ),
    paste0(
      "Completed at: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S %Z"
      )
    ),
    paste0(
      "QC: ",
      qc_file
    )
  )

  writeLines(
    completion_text,
    completion_file
  )

  cat(
    "\n============================================================\n",
    "VIRGINIA HPMS AADT 2011-2020 COMPLETED\n",
    "============================================================\n",
    sep = ""
  )

  print(
    qc[, .(
      year,
      source_type,
      n_features,
      total_field,
      combination_field,
      single_unit_field,
      percent_missing_total_aadt,
      total_aadt_median,
      total_aadt_p95,
      total_aadt_maximum,
      percent_missing_truck_aadt,
      truck_aadt_median,
      truck_aadt_p95,
      truck_aadt_maximum
    )]
  )

  cat(
    "\nOutput directory:\n",
    out_dir,
    "\nCompletion marker:\n",
    completion_file,
    "\n",
    sep = ""
  )
}

status <- tryCatch(
  {
    run_download()
    TRUE
  },
  error = function(e) {
    message(
      "\n============================================================\n",
      "HPMS AADT DOWNLOAD ABORTED\n",
      "============================================================\n",
      conditionMessage(e),
      "\n"
    )
    FALSE
  }
)

if (!isTRUE(status)) {
  if (interactive()) {
    stop(
      "HPMS AADT download stopped. Read the error above.",
      call. = FALSE
    )
  } else {
    quit(
      save = "no",
      status = 1L,
      runLast = FALSE
    )
  }
}
