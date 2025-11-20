# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# Incremental update (baseline 1990–2024 + rolling 2025–2029)

# 0 - Load metadata once (as in your working script) ----
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

# 1 - Helpers (parser + standardizer) ----
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

coerce_and_select_schema <- function(DT) {
  # Use the same column set/types every time to keep schemas aligned
  keep_cols <- c(
    "series_id","survey","seasonal_code","area_code","measure_code","area_type_code",
    "year","period","month","value",
    "area_text","series_title","measure_text","areatype_text","seasonal",
    "begin_year","begin_period","end_year","end_period","srd_code","srd_text"
  )
  keep_cols <- intersect(keep_cols, names(DT))
  DT <- DT[, ..keep_cols]
  
  char_cols <- intersect(c("series_id","survey","seasonal_code","area_code","measure_code",
                           "area_type_code","period",
                           "area_text","series_title","measure_text","areatype_text","seasonal",
                           "begin_period","end_period","srd_code","srd_text"),
                         names(DT))
  int_cols  <- intersect(c("year","month","begin_year","end_year"), names(DT))
  dbl_cols  <- intersect(c("value"), names(DT))
  
  for (cc in char_cols) DT[, (cc) := as.character(get(cc))]
  for (cc in int_cols)  DT[, (cc) := as.integer(get(cc))]
  for (cc in dbl_cols)  DT[, (cc) := as.numeric(get(cc))]
  
  DT
}

# 2 - Process the 2025–2029 raw file into "updates" parquet (monthly/annual) ----
process_updates_2529 <- function(
    raw_path_2529 = "data/BLSraw/bls_la.data.0.CurrentU25-29",
    out_dir_month_updates  = "data/BLSprocessed/monthly_updates_25_29",
    out_dir_annual_updates = "data/BLSprocessed/annual_updates_25_29"
) {
  DT <- utils::read.delim(raw_path_2529, sep = "\t", quote = "", stringsAsFactors = FALSE, header = TRUE)
  setDT(DT)
  
  DT[, series_id := gsub(" ", "", trimws(series_id))]
  parts <- parse_la_series_id(DT$series_id)
  DT[, c("survey","seasonal_code","area_code","measure_code") := parts]
  DT[, month := as.integer(sub("^M", "", period))]
  DT[period == "M13", month := NA_integer_]
  
  # --- metadata joins ---
  setkey(DT, measure_code);  DT <- measure[DT]
  setkey(DT, seasonal_code); DT <- seasonal[DT]
  setkey(DT, area_code);     DT <- area[, .(area_type_code, area_code, area_text)][DT]
  setkey(DT, area_type_code);DT <- area_type[DT]
  setkey(DT, period);        DT <- period_md[DT]
  setkey(DT, series_id);     DT <- series_md[, .(series_id, srd_code, series_title, begin_year, begin_period, end_year, end_period)][DT]
  setkey(DT, srd_code);      DT <- srd[DT]
  
  DT <- coerce_and_select_schema(DT)
  
  # split & write (overwrite the updates folders each run)
  DT_month  <- DT[period != "M13"]
  DT_annual <- DT[period == "M13"]
  
  unlink(out_dir_month_updates,  recursive = TRUE, force = TRUE)
  unlink(out_dir_annual_updates, recursive = TRUE, force = TRUE)
  dir.create(out_dir_month_updates,  recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_annual_updates, recursive = TRUE, showWarnings = FALSE)
  
  write_parquet(DT_month,  file.path(out_dir_month_updates,  paste0(basename(raw_path_2529), ".parquet")))
  write_parquet(DT_annual, file.path(out_dir_annual_updates, paste0(basename(raw_path_2529), ".parquet")))
}

# 3 - Combine baseline (<=2024) with updates (>=2025) into a NEW versioned output ----
combine_baseline_with_updates <- function(
    baseline_month_dir  = "data/BLSprocessed/monthly_parquet",
    baseline_annual_dir = "data/BLSprocessed/annual_parquet",
    updates_month_dir   = "data/BLSprocessed/monthly_updates_25_29",
    updates_annual_dir  = "data/BLSprocessed/annual_updates_25_29",
    out_month_dir_root  = "data/BLSprocessed/combined_monthly_parquet",
    out_annual_dir_root = "data/BLSprocessed/combined_annual_parquet"
) {
  # open datasets (unified schema)
  ds_base_m <- open_dataset(baseline_month_dir,  format = "parquet", unify_schemas = TRUE)
  ds_base_a <- open_dataset(baseline_annual_dir, format = "parquet", unify_schemas = TRUE)
  ds_up_m   <- open_dataset(updates_month_dir,   format = "parquet", unify_schemas = TRUE)
  ds_up_a   <- open_dataset(updates_annual_dir,  format = "parquet", unify_schemas = TRUE)
  
  # filter baseline to <=2024, updates to >=2025; union_all keeps latest 2025+ fully
  combo_m <- ds_base_m %>%
    filter(year <= 2024) %>%
    union_all(ds_up_m %>% filter(year >= 2025))
  
  combo_a <- ds_base_a %>%
    filter(year <= 2024) %>%
    union_all(ds_up_a %>% filter(year >= 2025))
  
  # versioned output dirs (e.g., by YYYYMMDD)
  stamp <- format(Sys.Date(), "%Y%m%d")
  out_month_dir  <- file.path(out_month_dir_root,  paste0("v", stamp))
  out_annual_dir <- file.path(out_annual_dir_root, paste0("v", stamp))
  
  # write partitioned datasets (fast; no giant in-memory collect)
  # overwrite version folder if re-run same day
  unlink(out_month_dir,  recursive = TRUE, force = TRUE)
  unlink(out_annual_dir, recursive = TRUE, force = TRUE)
  
  write_dataset(combo_m, out_month_dir,  format = "parquet", existing_data_behavior = "overwrite")
  write_dataset(combo_a, out_annual_dir, format = "parquet", existing_data_behavior = "overwrite")
  
  message("Wrote combined monthly to: ", out_month_dir)
  message("Wrote combined annual  to: ", out_annual_dir)
  
  invisible(list(month_dir = out_month_dir, annual_dir = out_annual_dir))
}


# 4 - Combine them all ----
  # Rebuild the updates parquet from the new 25-29 raw file
  process_updates_2529("data/BLSraw/bls_la.data.0.CurrentU25-29")

  # Combine baseline (<=2024) + updates (>=2025) into a NEW versioned folder
  paths <- combine_baseline_with_updates()

  
# 5 - Check exercise ----
  
  ds_month  <- open_dataset("data/BLSprocessed/combined_monthly_parquet", format = "parquet", unify_schemas = TRUE)
  
  seasonal_pref <- "U"  # "S" (adjusted), "U" (unadjusted), or NULL for both
  
  df_states <- ds_month %>%
    filter(
      measure_code == "03",
      area_type_code == "A",
      year >= 1990, year <= 2025
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
  
# 6 - Create final dataset monthly by county ----
  
  ds_month  <- open_dataset("data/BLSprocessed/combined_monthly_parquet", 
                            format = "parquet", unify_schemas = TRUE)
  
  ctowns <- read_excel("data/college_towns.xlsx",sheet = 2) %>% na.omit()
  
  
  df_counties <- ds_month %>%
    filter(
      measure_code == "03",
      area_type_code == "F",
      year >= 1990, year <= 2025
    ) %>%
    { if (is.null(seasonal_pref)) . else filter(., seasonal_code == seasonal_pref) } %>%
    mutate(
      date  = as.Date(paste0(year,"-", month, "-01"),format="%Y-%m-%d"),
      value = cast(value, float64(), safe = FALSE)
    ) %>%
    collect()
  
  df_counties <- df_counties %>% 
    mutate(
      state_fips  = as.integer(substr(area_code, 3, 4)),
      county_fips = as.integer(substr(area_code, 3, 7))
      ) %>% 
    mutate(county_fips = as.character(county_fips))
  
  ctowns2 <- ctowns %>%
    mutate(FIPS_county = as.character(FIPS_county)) %>%
    distinct(FIPS_county) %>% 
    mutate(CTind = 1)
  
  df_counties2 <- df_counties %>% 
    left_join(ctowns2, by = c("county_fips" = "FIPS_county")) %>%
    mutate(
      CTind = ifelse(is.na(CTind),0, 1)
    )
  
  # SA Data
  library(feasts)    # for STL decomposition with tsibble
  library(tsibble)
  library(imputeTS)   # for robust interpolation of missing values
  
  df_ts <- df_counties2 %>%
    mutate(
      date = yearmonth(date)  # from lubridate & tsibble
    ) %>%
    as_tsibble(
      key = county_fips,
      index = date
    )
  
  df_ts <- df_ts %>%
    group_by(county_fips) %>%
    mutate(
      value_interp = na_interpolation(value, option = "linear")
    ) %>%
    ungroup()
  
  df_stl <- df_ts %>%
    model(STL(value_interp ~ season(window = "periodic"))) %>%
    components() %>% 
    mutate(
      # seasonally adjusted = observed - seasonal
      value_sa = season_adjust
    ) %>%
    as_tibble()
  
  df_counties2 <- df_counties2 %>%
    mutate(date = yearmonth(as.Date(date))) %>%
    left_join(df_stl %>% select(county_fips, date, value_sa), 
              by = c("county_fips", "date"))
  
  df_counties2 <- df_counties2 %>%
    mutate(date = as.Date(date))
  
  
  ggplot(df_counties2, aes(x = date, y = value_sa, group = county_fips)) +
    # background: all counties in grey
    geom_line(color = "grey80", alpha = 0.4) +
    # highlight: college towns in black
    geom_line(
      data = df_counties2 %>% dplyr::filter(CTind == 1),
      color = "black",
      alpha = 0.9,
      linewidth  = 0.6
    ) +
    labs(
      x = "Date",
      y = "Unemployment rate",
      title = "Unemployment time series: all counties (grey) vs college towns (black)"
    ) +
    theme_minimal()
  
  
  # 1. Summarise by date and college_town
  df_summary <- df_counties2 %>%
    group_by(date, CTind) %>%
    summarise(
      p25 = quantile(value_sa, 0.25, na.rm = TRUE),
      p50 = quantile(value_sa, 0.50, na.rm = TRUE),
      p75 = quantile(value_sa, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      CTind = factor(
        CTind,
        levels = c(0, 1),
        labels = c("Non-college counties", "College-town counties")
      )
    )
  
  # 2. Plot: ribbon = IQR, line = median
  g1county <- ggplot(df_summary,
         aes(x = date)) +
    geom_ribbon(aes(ymin = p25, ymax = p75,
                    color = CTind,
                    fill  = CTind),
                alpha = 0.20,
                color = NA) +
    geom_line(aes(y = p50,linetype = CTind,color = CTind),
              size = 1) + 
    scale_linetype_manual(values = c("dotted","solid"),
                          labels=c('Non-college counties', 'College-town counties'))+
    labs(
      x = "Date",
      y = "Unemployment rate",
      color = "",
      fill  = "",linetype="",
      #title = "Unemployment in college vs non-college counties",
      subtitle = "Median (line) with 25th–75th percentiles shaded"
    ) +
    theme_bw(base_size = 16)+theme(legend.position=c(.2,.75))+
    labs(x="")+
    scale_x_date(date_breaks = "4 years",date_labels = "%Y")
  
  ggsave(filename = "output/g1county.png",g1county,scale=1,units = "cm",width = 20,height=15)
  
  
  # Fast DID
  
  recessdf <- read_excel("data/macro/USRECM.xlsx",sheet=2) %>% 
    mutate(date = as.Date(date,format="%Y-%m-%d"))
  
  df_counties2 <- df_counties2 %>% merge(recessdf,by="date",all.x = T) %>% 
    mutate(month_fe = factor(format(date, "%m")))

  library(fixest)
  did_mod <- feols(
    value_sa ~ CTind * USRECM | county_fips + date,
    cluster = ~ county_fips,
    data = df_counties2
  )
  
  summary(did_mod)
  
  
  # Event
  # Identify recession start months (from peak→trough indicator)
  df2 <- df_counties2 %>%
    mutate(
      # ensure Date class; your format is "%Y-%m-%d"
      date = as.Date(date, format = "%Y-%m-%d")
    )
  
  # Create a global time index t = 1, 2, 3, ... for each unique date
  time_index <- df2 %>%
    distinct(date) %>%
    arrange(date) %>%
    mutate(t = row_number())
  
  # Join back to your main data
  df2 <- df2 %>%
    left_join(time_index, by = "date") %>%
    arrange(county_fips, date)
  
  # Recession info by date (USRECM should be same for all counties each date)
  rec_info <- df2 %>%
    distinct(date, t, USRECM) %>%
    arrange(t)
  
  # A recession starts when USRECM goes from 0 -> 1
  rec_info <- rec_info %>%
    mutate(rec_start = ifelse(USRECM == 1 & dplyr::lag(USRECM, default = 0) == 0, 1, 0))
  
  # Vector of t indices where recessions start
  rec_start_t <- rec_info %>%
    filter(rec_start == 1) %>%
    pull(t)
  
  rec_start_t
  
  
  if (length(rec_start_t) == 0) {
    stop("No recession starts found (USRECM never goes 0 -> 1).")
  }
  
  # Function that gives distance (in periods) to nearest recession start
  closest_rel_time <- function(t, rec_starts) {
    # rec_starts is the vector rec_start_t
    diffs <- t - rec_starts
    diffs[which.min(abs(diffs))]
  }
  
  # Compute rel_time for each row
  df2$rel_time <- vapply(df2$t, closest_rel_time, integer(1L), 
                         rec_starts = rec_start_t)
  
  # Trim / bin the event window and pick a valid reference period
  min_k <- -24
  max_k <-  24
  
  df2 <- df2 %>%
    mutate(
      rel_group = case_when(
        rel_time < min_k ~ min_k,
        rel_time > max_k ~ max_k,
        TRUE ~ rel_time
      )
    )
  
  # Check what values you actually have
  sort(unique(df2$rel_group))
  
  # last pre-recession event period present (e.g. -1, -2, etc.)
  base_period <- max(df2$rel_group[df2$rel_group < 0], na.rm = TRUE)
  base_period
  
  # Run the event-study regression (fixest)
  event_est <- feols(
    value_sa ~ i(rel_group, CTind, ref = base_period) | county_fips + t,
    cluster = ~ county_fips,
    data = df2
  )
  
  summary(event_est)
  
  g2county_all <- ggiplot(event_est)+
    labs(
      #title = "Event Study: Effect of Recessions on Unemployment in College Towns",
      title="",subtitle = "Estimation by county",
      x = "Periods relative to recession start",
      y = "Effect on unemployment rate (percentage points)")+theme_bw(base_size = 16)
  
  ggsave(filename = "output/g2county.png",g2county_all,scale=1,units = "cm",width = 20,height=15)
 
  # Run the event-study regression (fixest) without COVID
  event_est <- feols(
    value_sa ~ i(rel_group, CTind, ref = base_period) | county_fips + t,
    cluster = ~ county_fips,
    data = df2 %>% filter(year<2020)
  )
  
  summary(event_est)
  
  g2county_pre <- ggiplot(event_est)+
    labs(
      #title = "Event Study: Effect of Recessions on Unemployment in College Towns",
      title="",subtitle = "Estimation by county",
      x = "Periods relative to recession start",
      y = "Effect on unemployment rate (percentage points)")+theme_bw(base_size = 16)
  
  ggsave(filename = "output/g2county_pre.png",g2county_pre,scale=1,units = "cm",width = 20,height=15)
  


# 7 - Create final dataset monthly by MSA ----
  
  df_msas <- ds_month %>%
    filter(
      measure_code == "03",
      area_type_code == "B",
      year >= 1990, year <= 2025
    ) %>%
    { if (is.null(seasonal_pref)) . else filter(., seasonal_code == seasonal_pref) } %>%
    mutate(
      date  = as.Date(paste0(year,"-", month, "-01"),format="%Y-%m-%d"),
      value = cast(value, float64(), safe = FALSE)
    ) %>%
    collect()
  
  df_msas <- df_msas %>% 
    mutate(
      state_fips  = as.integer(substr(area_code, 3, 4)),
      msa_fips = as.integer(substr(area_code, 5, 9))
    ) %>% 
    mutate(msa_fips = as.character(msa_fips))
  
  ctowns2 <- ctowns %>%
    mutate(FIPS_MSA = as.character(FIPS_MSA)) %>%
    distinct(FIPS_MSA) %>% 
    mutate(CTind = 1)
  
  df_msas2 <- df_msas %>% 
    left_join(ctowns2, by = c("msa_fips" = "FIPS_MSA")) %>%
    mutate(
      CTind = ifelse(is.na(CTind),0, 1)
    )
  
  df_ts <- df_msas2 %>%
    mutate(
      date = yearmonth(date)  # from lubridate & tsibble
    ) %>%
    as_tsibble(
      key = msa_fips,
      index = date
    )
  
  df_ts <- df_ts %>%
    group_by(msa_fips) %>%
    mutate(
      value_interp = na_interpolation(value, option = "linear")
    ) %>%
    ungroup()
  
  df_stl <- df_ts %>%
    model(STL(value_interp ~ season(window = "periodic"))) %>%
    components() %>% 
    mutate(
      # seasonally adjusted = observed - seasonal
      value_sa = season_adjust
    ) %>%
    as_tibble()
  
  df_msas2 <- df_msas2 %>%
    mutate(date = yearmonth(as.Date(date))) %>%
    left_join(df_stl %>% select(msa_fips, date, value_sa), 
              by = c("msa_fips", "date"))
  
  df_msas2 <- df_msas2 %>%
    mutate(date = as.Date(date))
  
  ggplot(df_msas2, aes(x = date, y = value_sa, group = msa_fips)) +
    # background: all msas in grey
    geom_line(color = "grey80", alpha = 0.4) +
    # highlight: college towns in black
    geom_line(
      data = df_msas2 %>% dplyr::filter(CTind == 1),
      color = "black",
      alpha = 0.9,
      linewidth  = 0.6
    ) +
    labs(
      x = "Date",
      y = "Unemployment rate",
      title = "Unemployment time series: all msas (grey) vs college towns (black)"
    ) +
    theme_minimal()
  
  
  # 1. Summarise by date and college_town
  df_summary <- df_msas2 %>%
    group_by(date, CTind) %>%
    summarise(
      p25 = quantile(value_sa, 0.25, na.rm = TRUE),
      p50 = quantile(value_sa, 0.50, na.rm = TRUE),
      p75 = quantile(value_sa, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      CTind = factor(
        CTind,
        levels = c(0, 1),
        labels = c("Non-college MSAs", "College-town MSAs")
      )
    )
  
  # 2. Plot: ribbon = IQR, line = median
  g1msa <- ggplot(df_summary,
                     aes(x = date)) +
    geom_ribbon(aes(ymin = p25, ymax = p75,
                    color = CTind,
                    fill  = CTind),
                alpha = 0.20,
                color = NA) +
    geom_line(aes(y = p50,linetype = CTind,color = CTind),
              size = 1) + 
    scale_linetype_manual(values = c("dotted","solid"),
                          labels=c('Non-college MSAs', 'College-town MSAs'))+
    labs(
      x = "Date",
      y = "Unemployment rate",
      color = "",
      fill  = "",linetype="",
      #title = "Unemployment in college vs non-college counties",
      subtitle = "Median (line) with 25th–75th percentiles shaded"
    ) +
    theme_bw(base_size = 16)+theme(legend.position=c(.2,.75))+
    labs(x="")+
    scale_x_date(date_breaks = "4 years",date_labels = "%Y")
  
  ggsave(filename = "output/g1msa.png",g1msa,scale=1,units = "cm",width = 20,height=15)
  
  
  
  # Fast DID
  
  recessdf <- read_excel("data/macro/USRECM.xlsx",sheet=2) %>% 
    mutate(date = as.Date(date,format="%Y-%m-%d"))
  
  df_msas2 <- df_msas2 %>% merge(recessdf,by="date",all.x = T) %>% 
    mutate(month_fe = factor(format(date, "%m")))
  
  library(fixest)
  did_mod <- feols(
    value_sa ~ CTind * USRECM | msa_fips + date,
    cluster = ~ msa_fips,
    data = df_msas2 %>% dplyr::filter(year<2020)
  )
  
  summary(did_mod)
  
  
  # Event
  # Identify recession start months (from peak→trough indicator)
  df2 <- df_msas2 %>%
    mutate(
      # ensure Date class; your format is "%Y-%m-%d"
      date = as.Date(date, format = "%Y-%m-%d")
    )
  
  # Create a global time index t = 1, 2, 3, ... for each unique date
  time_index <- df2 %>%
    distinct(date) %>%
    arrange(date) %>%
    mutate(t = row_number())
  
  # Join back to your main data
  df2 <- df2 %>%
    left_join(time_index, by = "date") %>%
    arrange(msa_fips, date)
  
  # Recession info by date (USRECM should be same for all counties each date)
  rec_info <- df2 %>%
    distinct(date, t, USRECM) %>%
    arrange(t)
  
  # A recession starts when USRECM goes from 0 -> 1
  rec_info <- rec_info %>%
    mutate(rec_start = ifelse(USRECM == 1 & dplyr::lag(USRECM, default = 0) == 0, 1, 0))
  
  # Vector of t indices where recessions start
  rec_start_t <- rec_info %>%
    filter(rec_start == 1) %>%
    pull(t)
  
  rec_start_t
  
  
  if (length(rec_start_t) == 0) {
    stop("No recession starts found (USRECM never goes 0 -> 1).")
  }
  
  # Function that gives distance (in periods) to nearest recession start
  closest_rel_time <- function(t, rec_starts) {
    # rec_starts is the vector rec_start_t
    diffs <- t - rec_starts
    diffs[which.min(abs(diffs))]
  }
  
  # Compute rel_time for each row
  df2$rel_time <- vapply(df2$t, closest_rel_time, integer(1L), 
                         rec_starts = rec_start_t)
  
  # Trim / bin the event window and pick a valid reference period
  min_k <- -24
  max_k <-  24
  
  df2 <- df2 %>%
    mutate(
      rel_group = case_when(
        rel_time < min_k ~ min_k,
        rel_time > max_k ~ max_k,
        TRUE ~ rel_time
      )
    )
  
  # Check what values you actually have
  sort(unique(df2$rel_group))
  
  # last pre-recession event period present (e.g. -1, -2, etc.)
  base_period <- max(df2$rel_group[df2$rel_group < 0], na.rm = TRUE)
  base_period
  
  
  # Run the event-study regression (fixest)
  event_est <- feols(
    value_sa ~ i(rel_group, CTind, ref = base_period) | msa_fips + t,
    cluster = ~ msa_fips,
    data = df2
  )
  
  summary(event_est)
  
  
  g2msa_all <- ggiplot(event_est)+
    labs(
      #title = "Event Study: Effect of Recessions on Unemployment in College Towns",
      title="",subtitle = "Estimation by MSAs",
      x = "Periods relative to recession start",
      y = "Effect on unemployment rate (percentage points)")+theme_bw(base_size = 16)
  
  ggsave(filename = "output/g2msa_all.png",g2msa_all,scale=1,units = "cm",width = 20,height=15)
  
  # Run the event-study regression (fixest)
  event_est <- feols(
    value_sa ~ i(rel_group, CTind, ref = base_period) | msa_fips + t,
    cluster = ~ msa_fips,
    data = df2 %>% dplyr::filter(year<2020)
  )
  
  summary(event_est)
  
  
  g2msa_pre <- ggiplot(event_est)+
    labs(
      #title = "Event Study: Effect of Recessions on Unemployment in College Towns",
      title="",subtitle = "Estimation by MSAs",
      x = "Periods relative to recession start",
      y = "Effect on unemployment rate (percentage points)")+theme_bw(base_size = 16)
  
  ggsave(filename = "output/g2msa_pre.png",g2msa_pre,scale=1,units = "cm",width = 20,height=15)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # --- cross-county distribution per state-month ---
  q_by_state <- df_counties2 %>%
    group_by(srd_text, date) %>%
    summarize(
      q05 = quantile(value, 0.05, na.rm = TRUE),
      q32 = quantile(value, 0.32, na.rm = TRUE),
      q50 = quantile(value, 0.50, na.rm = TRUE),
      q68 = quantile(value, 0.68, na.rm = TRUE),
      q95 = quantile(value, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # --- plot ribbons: 5–95 (lighter), 32–68 (darker), median line ---
  p <- q_by_state %>% filter(srd_text%in%c("Texas","New York","California","Florida",
                                           "Illinois","Pennsylvania","Utah")) %>% ggplot( aes(x = date)) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.20) +
    geom_ribbon(aes(ymin = q32, ymax = q68), alpha = 0.35) +
    geom_line(aes(y = q50), linewidth = 0.35) +
    facet_wrap(~ srd_text, ncol = 3, scales = "free_y") +
    scale_x_date(labels = date_format("%Y"), breaks = pretty_breaks(6)) +
    labs(
      title = "Distribution of County Unemployment Rates by State",
      subtitle = sprintf("%d–%d (%s)", 1990, 2025,
                         if (is.null(seasonal_pref)) "All seasonal statuses"
                         else if (seasonal_pref == "S") "Seasonally Adjusted"
                         else "Not Seasonally Adjusted"),
      x = NULL, y = "Percent",
      caption = "Source: BLS LAUS; ribbons show 5–95 and 32–68 quantile bands; line is median."
    ) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8),
          plot.title = element_text(face = "bold"))
  
  print(p)
  