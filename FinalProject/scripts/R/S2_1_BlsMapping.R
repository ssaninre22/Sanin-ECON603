# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# 0 - Helpers ----
# Minimal, dependency-free reader for tab-delimited BLS mapping tables.
# - Keeps text as character (not factors)
# - Handles UTF-8 / UTF-8-BOM / latin1 gracefully
# - Trims stray whitespace from both columns
read_map <- function(path,expected_name = NULL,col_class = c("character", "character")) {
  try_read <- function(enc) {
    utils::read.delim(
      file = path,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      quote = "",            # do not treat quotes specially
      check.names = FALSE,   # keep original headers if needed
      fileEncoding = enc,
      colClasses = col_class
    )
  }
  
  df <- tryCatch(try_read("UTF-8"),
                 error = function(e) tryCatch(try_read("UTF-8-BOM"),
                                              error = function(e) try_read("latin1")))
  
  # Trim whitespace from character columns
  for (j in seq_along(df)) {
    if (is.character(df[[j]])) df[[j]] <- trimws(df[[j]])
  }
  
  # Optionally enforce expected column names
  if (!is.null(expected_name)) {
    if (length(expected_name) != ncol(df)) {
      stop(sprintf("Expected %d columns but found %d in %s",
                   length(expected_name), ncol(df), path))
    }
    names(df) <- expected_name
  }
  
  df
}

# 1 - Read mapping tables ----

  # Paths
  path_area_type <- "data/BLSmeta/la.area_type.txt"
  path_area      <- "data/BLSmeta/la.area.txt"
  path_measure   <- "data/BLSmeta/la.measure.txt"
  path_period    <- "data/BLSmeta/la.period"
  path_seasonal  <- "data/BLSmeta/la.seasonal.txt"
  path_footnote  <- "data/BLSmeta/la.footnote.txt"
  path_series    <- "data/BLSmeta/la.series.txt"
  path_srd       <- "data/BLSmeta/la.state_region_division.txt"
  
  # Colnames
  cname_area_type <- c("area_type_code", "areatype_text")
  cname_area      <- c("area_type_code", "area_code","area_text",	"display_level",
                       "selectable", "sort_sequence")
  cname_measure   <- c("measure_code", "measure_text")
  cname_period    <- c("period", "period_abbr","period_name")
  cname_seasonal  <- c("seasonal_code", "seasonal_text")
  cname_footnote  <- c("footnote_code", "footnote_text")
  cname_series    <- c("series_id", "area_type_code",	"area_code","measure_code",	
                       "seasonal",	"srd_code",	"series_title",	"footnote_codes",	
                       "begin_year",	"begin_period",	"end_year",	"end_period")
  cname_srd       <- c("srd_code", "srd_text")
  
  # Coltypes
  ctype_area_type <- c(rep("character",2))
  ctype_area      <- c(rep("character",5),"numeric") 
  ctype_measure   <- c(rep("character",2))
  ctype_period    <- c(rep("character",3))
  ctype_seasonal  <- c(rep("character",2))
  ctype_footnote  <- c(rep("character",2))
  ctype_series    <- c(rep("character",8),"numeric","character","numeric","character")
  ctype_srd       <- c(rep("character",2))
  
  
  # The BLS header row in your file is:
  area_type <- read_map(path_area_type,expected_name = cname_area_type, col_class = ctype_area_type)
  area      <- read_map(path_area     ,expected_name = cname_area     , col_class = ctype_area     )
  measure   <- read_map(path_measure  ,expected_name = cname_measure  , col_class = ctype_measure  )
  period    <- read_map(path_period   ,expected_name = cname_period   , col_class = ctype_period   )
  seasonal  <- read_map(path_seasonal ,expected_name = cname_seasonal , col_class = ctype_seasonal )
  footnote  <- read_map(path_footnote ,expected_name = cname_footnote , col_class = ctype_footnote )
  series    <- read_map(path_series   ,expected_name = cname_series   , col_class = ctype_series   )
  srd       <- read_map(path_srd      ,expected_name = cname_srd      , col_class = ctype_srd      )
  
  
  # Save as RData (.rda)
  saveRDS( area_type , "data/BLSmeta/la_area_type.rds" )
  saveRDS( area      , "data/BLSmeta/la_area.rds" )
  saveRDS( measure   , "data/BLSmeta/la_measure.rds")
  saveRDS( period    , "data/BLSmeta/la_period.rds")
  saveRDS( seasonal  , "data/BLSmeta/la_seasonal.rds")
  saveRDS( footnote  , "data/BLSmeta/la_footnote.rds")
  saveRDS( series    , "data/BLSmeta/la_series.rds")      
  saveRDS( srd       , "data/BLSmeta/la_srd.rds")


