# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# 0 - Load all the needed datasets ----

  # 0.1 Load all the metada files
  area_type <- readRDS("data/BLSmeta/la_area_type.rds")
  area      <- readRDS("data/BLSmeta/la_area.rds" )
  measure   <- readRDS("data/BLSmeta/la_measure.rds")
  period    <- readRDS("data/BLSmeta/la_period.rds")
  seasonal  <- readRDS("data/BLSmeta/la_seasonal.rds")
  footnote  <- readRDS("data/BLSmeta/la_footnote.rds")
  series    <- readRDS("data/BLSmeta/la_series.rds")      
  srd       <- readRDS("data/BLSmeta/la_srd.rds")
  
  # 0.2 - Combine all the old data sets 90-24
  df1 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU90-94")
  df2 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU95-99")
  df3 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU00-04")
  df4 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU05-09")
  df5 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU10-14")
  df6 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU15-19")
  df7 <- read.delim("data/BLSraw/bls_la.data.0.CurrentU20-24")
  
  dfold <- rbind(df1,df2,df3,df4,df5,df6,df7)
  saveRDS(dfold,"data/BLSprocessed/Current9024.rds")
  
  rm(df1,df2,df3,df4,df5,df6,df7,dfold)
  gc()
  
  # 0.3 - Combine Old file with metadata
    # Utility function:
      # Parse LAS series_id into components
      parse_la_series_id <- function(x) {
        x <- trimws(x)
        if (any(nchar(x) != 20L)) stop("All series_id must have length 20.")
        data.frame(
          survey        = substr(x, 1, 2),
          seasonal_code = substr(x, 3, 3),
          area_code     = substr(x, 4, 18),
          measure_code  = substr(x, 19, 20),
          stringsAsFactors = FALSE
        )
      }
      
      combine_df_meta <- function(namefile){
        df <- read.delim(paste0("data/BLSraw/bls_la.data.0.",namefile))
        # Parse
        df_part <- parse_la_series_id(df$series_id)
        df <- cbind(df, df_part)
        df$series_id <-  gsub(" ", "", df$series_id)
        # Merger
        df_temp <- df %>% 
          merge(measure, by = "measure_code", all.x = TRUE) %>% 
          merge(seasonal, by = "seasonal_code", all.x = TRUE) %>% 
          merge(area %>% select(area_type_code,area_code,area_text),by = "area_code",all.x=TRUE) %>% 
          merge(area_type, by = "area_type_code",all.x = TRUE) %>% 
          merge(period, by = "period",all.x=TRUE) %>% 
          merge(series %>% 
                  select(series_id,srd_code,series_title,begin_year,begin_period,end_year,end_period),
                by = "series_id",all.x = TRUE) %>% 
          merge(srd,by = "srd_code",all.x = TRUE)
        
        # Define month and annual datasets
        df_temp_mth <- df_temp %>% filter(period!="M13") %>% 
          mutate(month = as.character(str_remove(period, "M"))) %>% 
          mutate(date = as.Date(paste0(year,"-",month,"-01"),format="%Y-%m-%d"))%>% 
          mutate(value = as.numeric(value))
        
        df_temp_annual <- df_temp %>% filter(period=="M13") %>% 
          mutate(date = as.Date(paste0(year,"-",12,"-01"),format="%Y-%m-%d"))%>% 
          mutate(value = as.numeric(value))
        
        saveRDS(df_temp_mth,paste0("data/BLSprocessed/",namefile,"_m.rds"))
        saveRDS(df_temp_annual,paste0("data/BLSprocessed/",namefile,"_y.rds"))
      }
      
      # Run and save
      combine_df_meta( "CurrentU90-94" ) # Ready
      combine_df_meta( "CurrentU95-99" ) # Ready
      combine_df_meta( "CurrentU00-04" ) # Ready
      combine_df_meta( "CurrentU05-09" ) # 
      combine_df_meta( "CurrentU10-14" )
      combine_df_meta( "CurrentU15-19" )
      combine_df_meta( "CurrentU20-24" )
      
     cs_data <- df_temp2 %>% filter(grepl("CT4815976000000", series_id, ignore.case = TRUE))



