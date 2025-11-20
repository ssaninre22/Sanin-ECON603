# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# 0 - Base folders ----
  qcew_raw_dir  <- "data/QCEWraw"
  qcew_proc_dir <- "data/QCEWprocessed"
  
  if (!dir.exists(qcew_proc_dir)) {
    dir.create(qcew_proc_dir, recursive = TRUE)
  }
  
  ## Standard column names for 1990-2025 files ----
  qcew_cols_1990on <- c(
    "area_fips","own_code","industry_code","agglvl_code","size_code","year","qtr",
    "disclosure_code","area_title","own_title","industry_title","agglvl_title",
    "size_title","qtrly_estabs_count","month1_emplvl","month2_emplvl","month3_emplvl",
    "total_qtrly_wages","taxable_qtrly_wages","qtrly_contributions","avg_wkly_wage",
    "lq_disclosure_code","lq_qtrly_estabs_count","lq_month1_emplvl","lq_month2_emplvl",
    "lq_month3_emplvl","lq_total_qtrly_wages","lq_taxable_qtrly_wages",
    "lq_qtrly_contributions","lq_avg_wkly_wage","oty_disclosure_code",
    "oty_qtrly_estabs_count_chg","oty_qtrly_estabs_count_pct_chg",
    "oty_month1_emplvl_chg","oty_month1_emplvl_pct","oty_month2_emplvl_chg",
    "oty_month2_emplvl_pct","oty_month3_emplvl_chg","oty_month3_emplvl_pct",
    "oty_total_qtrly_wages_chg","oty_total_qtrly_wages_pct",
    "oty_taxable_qtrly_wages_chg","oty_taxable_qtrly_wages_pct",
    "oty_qtrly_contributions_chg","oty_qtrly_contributions_pct",
    "oty_avg_wkly_wage_chg","oty_avg_wkly_wage_pct"
  )


# 1 - Processing raw data files ----
  # 1.0 - Helpers ----
    read_qcew_file <- function(path, year) {
      
      if (year <= 1989) {
        # 1975-1989: SIC, no locality/delta columns, names are already fine
        df <- readr::read_csv(
          file      = path,
          col_types = cols(.default = col_character())
        )
        
      } else if (year >= 1990 && year <= 2015) {
        # 1990-2015: NAICS, duplicated column name in header
        # -> skip header row and assign correct names
        df <- readr::read_csv(
          file      = path,
          col_types = cols(.default = col_character()),
          col_names = qcew_cols_1990on,
          skip      = 1
        )
        
      } else {
        # 2016-2025: NAICS, *_pct_chg names that must become *_pct
        df <- readr::read_csv(
          file      = path,
          col_types = cols(.default = col_character())
        )
        
        # Rename *_pct_chg -> *_pct to match 1990-2015 convention
        new_names <- names(df) %>%
          stringr::str_replace("_pct_chg$", "_pct")
        
        # Ensure they line up with qcew_cols_1990on (for 47-col files)
        if (length(new_names) == length(qcew_cols_1990on)) {
          names(df) <- new_names
        } else {
          warning("Unexpected number of columns in file: ", path,
                  " (", length(new_names), " columns)")
          names(df) <- new_names
        }
      }
      
      # Keep the source file name for debugging / tracing if needed
      df$source_file <- basename(path)
      df
    }
  
  # 1.1 - Process one QCEW year ----
    process_qcew_year <- function(year,
                                  raw_root = qcew_raw_dir) {
      
      message("Processing year: ", year)
      
      year_dir <- file.path(raw_root, year)
      
      if (!dir.exists(year_dir)) {
        warning("Directory for year ", year, " not found: ", year_dir)
        return(invisible(NULL))
      }
      
      csv_files <- list.files(
        path       = year_dir,
        pattern    = "\\.csv$",
        full.names = TRUE
      )
      
      if (length(csv_files) == 0) {
        warning("No CSV files found for year ", year, " in ", year_dir)
        return(invisible(NULL))
      }
      
      # Read and row-bind all industry/subindustry files for this year
      year_df <- purrr::map_dfr(
        .x = csv_files,
        .f = ~ read_qcew_file(.x, year = year)
      )
      
      # Save as RDS in the year folder
      out_path <- file.path(year_dir, paste0("CombinedIndustries", year, ".rds"))
      saveRDS(year_df, out_path)
      message("Saved: ", out_path,
              "  (", nrow(year_df), " rows, ", ncol(year_df), " cols)")
      
      # Clean memory
      rm(year_df)
      gc()
      
      invisible(out_path)
    }
  
  # 1.2 - Combine yearly Combined_industries ----
    
    # Main combination function
      combine_qcew_years <- function(years,
                                     out_filename,
                                     raw_root  = qcew_raw_dir,
                                     proc_root = qcew_proc_dir) {
        
        stopifnot(length(years) > 0)
        
        message("Combining years: ", paste(range(years), collapse = "-"),
                " into ", out_filename)
        
        rds_paths <- file.path(
          raw_root,
          years,
          paste0("CombinedIndustries", years, ".rds")
        )
        
        # Check which exist
        exists_vec <- file.exists(rds_paths)
        if (!all(exists_vec)) {
          warning("Some CombinedIndustries{year}.rds files are missing:\n",
                  paste(rds_paths[!exists_vec], collapse = "\n"))
        }
        
        rds_paths <- rds_paths[exists_vec]
        if (length(rds_paths) == 0) {
          warning("No RDS files found for years: ", paste(years, collapse = ", "))
          return(invisible(NULL))
        }
        
        # Read and bind
        df_list <- lapply(rds_paths, readRDS)
        combined_df <- dplyr::bind_rows(df_list)
        
        out_path <- file.path(proc_root, out_filename)
        saveRDS(combined_df, out_path)
        
        message("Saved combined file: ", out_path,
                "  (", nrow(combined_df), " rows, ", ncol(combined_df), " cols)")
        
        rm(df_list, combined_df)
        gc()
        
        invisible(out_path)
      }
    
    # Wrapper:
    build_all_qcew_groups <- function() {
      # 1975–1989
      combine_qcew_years(years = 1975:1989,
                         out_filename = "CombinedIndustries1975to1989.rds")
      
      # 1990–1999
      combine_qcew_years(years = 1990:1999,
                         out_filename = "CombinedIndustries1990to1999.rds")
      
      # 2000–2009
      combine_qcew_years(years = 2000:2009,
                         out_filename = "CombinedIndustries2000to2009.rds")
      
      # 2010–2019
      combine_qcew_years(years = 2010:2019,
                         out_filename = "CombinedIndustries2010to2019.rds")
      
      # 2020–2025
      combine_qcew_years(years = 2020:2025,
                         out_filename = "CombinedIndustries2020to2025.rds")
    }
    
  # 1.3 - Execute raw data readers ----
    # Step 1: create CombinedIndustries{yyyy}.rds inside each year folder
      years_to_process <- 2016:2025
      purrr::walk(years_to_process, process_qcew_year)
      
    # Step 2: create decade / period panels in data/QCEWprocessed/
      build_all_qcew_groups()
      
      
      # TESTS:
      
      xx <- readRDS(paste0(qcew_proc_dir,"/CombinedIndustries1975to1989.rds"))
      
      xy <- xx %>% filter(agglvl_code == "26") %>% 
        dplyr::select(area_fips,industry_code,year,qtr,area_title,industry_title,month1_emplvl,
                      month2_emplvl,month3_emplvl)
        
      
      xy %>% filter(grepl("Brazos", area_title)) %>% 
        mutate(year = as.numeric(year),
               qtr= as.numeric(qtr),
               month1_emplvl = as.numeric(month1_emplvl),
               month2_emplvl = as.numeric(month2_emplvl),
               month3_emplvl = as.numeric(month3_emplvl)) %>% 
        mutate(date = as.Date(paste0(year,"-",qtr*3,"-",1),format="%Y-%m-%d")) %>% 
        ggplot(aes(x=date,y=month3_emplvl))+geom_line()+theme_bw()
      
      
      xx <- readRDS(paste0(qcew_proc_dir,"/CombinedIndustries1990to1999.rds"))
      
      xx %>% dplyr::select(agglvl_code,agglvl_title) %>% unique() %>% print(n=30)
      
      xy <- xx %>% filter(agglvl_code == "70") %>% 
        dplyr::select(area_fips,industry_code,year,qtr,area_title,industry_title,month1_emplvl,
                      month2_emplvl,month3_emplvl)
    
      xy %>% filter(grepl("Santa Cruz", area_title)) %>% 
        mutate(year = as.numeric(year),
               qtr= as.numeric(qtr),
               month1_emplvl = as.numeric(month1_emplvl),
               month2_emplvl = as.numeric(month2_emplvl),
               month3_emplvl = as.numeric(month3_emplvl)) %>% 
        mutate(date = as.Date(paste0(year,"-",qtr*3,"-",1),format="%Y-%m-%d")) %>% 
        ggplot(aes(x=date,y=month3_emplvl))+geom_line()+theme_bw()
      
      
      xx <- readRDS(paste0(qcew_proc_dir,"/CombinedIndustries2000TO2009.rds"))
      
      xy <- xx %>% filter(industry_code=="61",grepl("Santa Cruz", area_title))
      
      
      xy %>% filter(own_code=="3") %>% 
        mutate(year = as.numeric(year),
               qtr= as.numeric(qtr),
               month1_emplvl = as.numeric(month1_emplvl),
               month2_emplvl = as.numeric(month2_emplvl),
               month3_emplvl = as.numeric(month3_emplvl)) %>% 
        mutate(date = as.Date(paste0(year,"-",qtr*3,"-",1),format="%Y-%m-%d")) %>% 
        ggplot(aes(x=date,y=month3_emplvl))+geom_line()+theme_bw()
      
      xx <- readRDS(paste0(qcew_proc_dir,"/CombinedIndustries2010TO2019.rds"))
      
      xy <- xx %>% filter(agglvl_code == "70") %>% 
        dplyr::select(area_fips,industry_code,year,qtr,area_title,industry_title,month1_emplvl,
                      month2_emplvl,month3_emplvl)
      
      
      xy %>% filter(grepl("Brazos", area_title)) %>% 
        mutate(year = as.numeric(year),
               qtr= as.numeric(qtr),
               month1_emplvl = as.numeric(month1_emplvl),
               month2_emplvl = as.numeric(month2_emplvl),
               month3_emplvl = as.numeric(month3_emplvl)) %>% 
        mutate(date = as.Date(paste0(year,"-",qtr*3,"-",1),format="%Y-%m-%d")) %>% 
        ggplot(aes(x=date,y=month3_emplvl))+geom_line()+theme_bw()
      