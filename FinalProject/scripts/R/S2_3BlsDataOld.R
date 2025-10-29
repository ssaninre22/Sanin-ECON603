# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# 0 - Load metadata and keys ----
  area_type <- readRDS("data/BLSmeta/la_area_type.rds") |> as.data.table()
  area      <- readRDS("data/BLSmeta/la_area.rds"     ) |> as.data.table()
  measure   <- readRDS("data/BLSmeta/la_measure.rds"  ) |> as.data.table()
  period_md <- readRDS("data/BLSmeta/la_period.rds"   ) |> as.data.table()
  seasonal  <- readRDS("data/BLSmeta/la_seasonal.rds" ) |> as.data.table()
  series_md <- readRDS("data/BLSmeta/la_series.rds"   ) |> as.data.table()
  srd       <- readRDS("data/BLSmeta/la_srd.rds"      ) |> as.data.table()
  
  setkey(area_type, area_type_code)
  setkey(area,      area_code)
  setkey(measure,   measure_code)
  setkey(period_md, period)
  setkey(seasonal,  seasonal_code)
  setkey(series_md, series_id)
  setkey(srd,       srd_code)

# 1 - Functions: utilities and main functions ----

  parse_la_series_id <- function(x) {
    x <- trimws(x)
    if (any(nchar(x) != 20L)) stop("All series_id must have length 20.")
    data.table(
      survey        = substr(x, 1, 2),
      seasonal_code = substr(x, 3, 3),
      area_code     = substr(x, 4, 18),
      measure_code  = substr(x, 19, 20)
    )
  }

  # 1.1 - Chunk processor: read -> parse -> join -> coerce -> write
  process_chunk <- function(path_in,
                            out_dir_month  = "data/BLSprocessed/monthly_parquet",
                            out_dir_annual = "data/BLSprocessed/annual_parquet") {
    
    # Read raw tab-delimited (keep spaces within fields)
    DT <- utils::read.delim(path_in, sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE)
    setDT(DT)
    
    # Normalize/parse identifiers
    DT[, series_id := gsub(" ", "", trimws(series_id))]
    parts <- parse_la_series_id(DT$series_id)
    DT[, c("survey","seasonal_code","area_code","measure_code") := parts]
    
    # Month from "M01".."M12"; keep NA for "M13" (annual)
    DT[, month := as.integer(sub("^M", "", period))]
    DT[period == "M13", month := NA_integer_]
    
    # ------------------
    # Metadata joins
    # ------------------
    # Join measure (adds measure_text, etc.)
    setkey(DT, measure_code)
    DT <- measure[DT]
    
    # Join seasonal (labels for S/U)
    setkey(DT, seasonal_code)
    DT <- seasonal[DT]
    
    # Join area (restrict to needed cols to avoid collisions)
    setkey(DT, area_code)
    DT <- area[, .(area_type_code, area_code, area_text)][DT]
    
    # Join area_type (human-readable area type)
    setkey(DT, area_type_code)
    DT <- area_type[DT]
    
    # Join period metadata (e.g., month names if available)
    setkey(DT, period)
    DT <- period_md[DT]
    
    # Join series metadata (title, begin/end, srd_code)
    setkey(DT, series_id)
    DT <- series_md[, .(series_id, srd_code, series_title, begin_year, begin_period, end_year, end_period)][DT]
    
    # Join srd metadata
    setkey(DT, srd_code)
    DT <- srd[DT]
    
    # ------------------
    # Coerce types BEFORE writing (fix schema mismatches)
    # ------------------
    # Keep one consistent column set and types
    keep_cols <- c(
      # identifiers / keys
      "series_id","survey","seasonal_code","area_code","measure_code","area_type_code",
      # time/value
      "year","period","month","value",
      # human-readable fields
      "area_text","series_title","measure_text","areatype_text","seasonal",
      # period/series metadata (if present)
      "begin_year","begin_period","end_year","end_period","srd_code","srd_text"
    )
    keep_cols <- intersect(keep_cols, names(DT))  # in case some meta fields don't exist
    
    DT <- DT[, ..keep_cols]
    
    # Coerce to stable schema
    char_cols <- intersect(c("series_id","survey","seasonal_code","area_code","measure_code",
                             "area_type_code","period",
                             "area_text","series_title","measure_text","areatype_text","seasonal",
                             "begin_period","end_period","srd_code","srd_text"),
                           names(DT))
    num_int_cols <- intersect(c("year","month","begin_year","end_year"), names(DT))
    num_dbl_cols <- intersect(c("value"), names(DT))
    
    for (cc in char_cols)     DT[, (cc) := as.character(get(cc))]
    for (cc in num_int_cols)  DT[, (cc) := as.integer(get(cc))]
    for (cc in num_dbl_cols)  DT[, (cc) := as.numeric(get(cc))]
    
    # Split monthly vs annual
    DT_month  <- DT[period != "M13"]
    DT_annual <- DT[period == "M13"]
    
    # Output dirs
    dir.create(out_dir_month,  recursive = TRUE, showWarnings = FALSE)
    dir.create(out_dir_annual, recursive = TRUE, showWarnings = FALSE)
    
    # Write parquet shards (consistent schema across all files)
    write_parquet(DT_month,  file.path(out_dir_month,  paste0(basename(path_in), ".parquet")))
    write_parquet(DT_annual, file.path(out_dir_annual, paste0(basename(path_in), ".parquet")))
    
    invisible(NULL)
  }

# 2 - Build monthly/annual parquet from raw chunks ----
  files <- c("data/BLSraw/bls_la.data.0.CurrentU90-94",
             "data/BLSraw/bls_la.data.0.CurrentU95-99",
             "data/BLSraw/bls_la.data.0.CurrentU00-04",
             "data/BLSraw/bls_la.data.0.CurrentU05-09",
             "data/BLSraw/bls_la.data.0.CurrentU10-14",
             "data/BLSraw/bls_la.data.0.CurrentU15-19",
             "data/BLSraw/bls_la.data.0.CurrentU20-24")

  # # Optional: start fresh by clearing output dirs
  # unlink("data/BLSprocessed/monthly_parquet", recursive = TRUE, force = TRUE)
  # unlink("data/BLSprocessed/annual_parquet",  recursive = TRUE, force = TRUE)

  invisible(lapply(files, process_chunk))
  gc()

# 2.1 Open combined datasets to check

  ds_month  <- open_dataset("data/BLSprocessed/monthly_parquet", format = "parquet", unify_schemas = TRUE)
  ds_annual <- open_dataset("data/BLSprocessed/annual_parquet",  format = "parquet", unify_schemas = TRUE)

  # Save into memory
   all_month  <- ds_month  %>% collect() %>% as.data.table()
   all_annual <- ds_annual %>% collect() %>% as.data.table()
  
   saveRDS(all_month,"data/BLSprocessed/old_month.rds")
   saveRDS(all_annual,"data/BLSprocessed/old_annual.rds")

# 2.2 Plot exercise: Unemployment rate by state, 1990–2024
  seasonal_pref <- "U"  # "S" (adjusted), "U" (unadjusted), or NULL for both
  
  df_states <- ds_month %>%
    filter(
      measure_code == "03",
      area_type_code == "A",
      year >= 1990, year <= 2024
    ) %>%
    { if (is.null(seasonal_pref)) . else filter(., seasonal_code == seasonal_pref) } %>%
    mutate(
      date  = as.Date(paste0(year,"-", month, "-01"),format="%Y-%m-%d"),
      value = cast(value, float64(), safe = FALSE)
    ) %>%
    select(area_text, area_code, date, value) %>%
    collect()
  
  # Plot (faceted by state)
  p <- ggplot(df_states, aes(date, as.numeric(value), group = area_code)) +
    geom_line(linewidth = 0.3) +
    facet_wrap(~ area_text, ncol = 6, scales = "free_y") +
    scale_x_date(labels = date_format("%Y"), breaks = pretty_breaks(6)) +
    labs(
      title = "Unemployment Rate by State (1990–2024)",
      subtitle = if (is.null(seasonal_pref)) "Seasonal: All"
      else if (seasonal_pref == "S") "Seasonally Adjusted"
      else "Not Seasonally Adjusted",
      x = NULL, y = "Percent",
      caption = "Source: BLS LAUS"
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8),
          plot.title = element_text(face = "bold"))
  
  print(p)




cs_data <- df_temp2 %>% filter(grepl("CT4815976000000", series_id, ignore.case = TRUE))
