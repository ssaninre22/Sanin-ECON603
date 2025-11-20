# ---------------------------------------------------------------------------- #
# Title: Universities as Macro-Stabilizers                                     #
# Author: Sebastian Sanin-Restrepo                                             #
# Code description: This is the main script that call all functions to obtain  #
#                   the output figures and tables in the document.             #
# ---------------------------------------------------------------------------- #

# install.packages(c("httr2", "fs"))
library(httr2)
library(fs)

BLS_LA_BASE <- "https://download.bls.gov/pub/time.series/la/"
BLS_TARGET_DIR <- "data/BLSraw"
BLS_MANIFEST <- file.path(BLS_TARGET_DIR, "_bls_manifest.csv")

# ---- 1) List remote files on the LAUS directory (public index) ----
# pattern example: "^la\\.data\\." for all la.data.* files
bls_la_list_remote <- function(pattern = NULL) {
  dir_url <- BLS_LA_BASE
  resp <- request(dir_url) |>
    req_user_agent("R/httr2 LAUS downloader (contact: your_email@example.com)") |>
    req_perform()
  
  html <- resp_body_string(resp)
  
  # Extract file names from hrefs in the Apache-style index
  links <- regmatches(
    html,
    gregexpr('href="([^"]+)"', html, perl = TRUE)
  )[[1]]
  files <- sub('^href="([^"]+)".*$', "\\1", links)
  
  # Keep only entries that look like files (skip parent dirs, icons, etc.)
  # The LA directory contains files like la.area, la.data.10.Arkansas, la.series, etc.
  files <- files[!grepl("^\\?|^/$|\\.\\./|^icons/|^$", files)]
  if (!is.null(pattern)) {
    files <- files[grepl(pattern, files)]
  }
  unname(files)
}

# ---- 2) Download a remote file if it is newer than our local copy ----
# Saves as data/raw/bls_{orig_name}
# Tracks ETag/Last-Modified in a small manifest CSV for smart updates
bls_la_download_if_new <- function(filename="la.data.19.Idaho",
                                   base_url = BLS_LA_BASE,
                                   target_dir = BLS_TARGET_DIR,
                                   prefix = "bls_",
                                   manifest_path = BLS_MANIFEST,
                                   retries = 3,
                                   verbose = TRUE) {
  dir_create(target_dir, recurse = TRUE)
  
  url <- paste0(base_url, filename)
  local_path <- file.path(target_dir, paste0(prefix, filename))
  
  # Load existing manifest (if any)
  man <- if (file_exists(manifest_path)) {
    tryCatch(read.csv(manifest_path, stringsAsFactors = FALSE),
             error = function(e) data.frame())
  } else data.frame()
  
  if (nrow(man) == 0 || !all(c("file","etag","last_modified","size_bytes","downloaded_at") %in% names(man))) {
    man <- data.frame(
      file = character(),
      etag = character(),
      last_modified = character(),
      size_bytes = numeric(),
      downloaded_at = character(),
      stringsAsFactors = FALSE
    )
  }
  
  # Prior metadata (if any)
  row_idx <- which(man$file == filename)
  prev_etag <- if (length(row_idx)) man$etag[row_idx] else NA_character_
  prev_lastmod <- if (length(row_idx)) man$last_modified[row_idx] else NA_character_
  
  # HEAD request to check freshness
  head_req <- request(url) |>
    req_user_agent("R/httr2 LAUS downloader (contact: your_email@example.com)") |>
    req_method("HEAD") |>
    req_retry(max_tries = retries)
  
  head_resp <- tryCatch(req_perform(head_req), error = function(e) e)
  if (inherits(head_resp, "error")) {
    if (verbose) message(sprintf("[WARN] HEAD failed for %s: %s", filename, head_resp$message))
    # Fallback: try to GET anyway (some servers block HEAD)
  }
  
  remote_etag <- if (!inherits(head_resp, "error")) resp_header(head_resp, "etag") else NA_character_
  remote_lastmod <- if (!inherits(head_resp, "error")) resp_header(head_resp, "last-modified") else NA_character_
  
  # If ETag or Last-Modified match our manifest, skip downloading
  if (!is.na(prev_etag) && !is.na(remote_etag) && identical(prev_etag, remote_etag)) {
    if (verbose) message(sprintf("[SKIP] %s unchanged (ETag match).", filename))
    return(invisible(local_path))
  }
  if (!is.na(prev_lastmod) && !is.na(remote_lastmod) && identical(prev_lastmod, remote_lastmod)) {
    if (verbose) message(sprintf("[SKIP] %s unchanged (Last-Modified match).", filename))
    return(invisible(local_path))
  }
  
  # Conditional GET with If-None-Match / If-Modified-Since (best-effort)
  get_req <- request(url) |>
    req_user_agent("R/httr2 LAUS downloader (contact: your_email@example.com)") |>
    req_retry(max_tries = retries)
  
  if (!is.na(prev_etag))       get_req <- req_headers(get_req, "If-None-Match" = prev_etag)
  if (!is.na(prev_lastmod))    get_req <- req_headers(get_req, "If-Modified-Since" = prev_lastmod)
  
  tmp <- tempfile("bls_", fileext = ".tmp")
  get_resp <- tryCatch(req_perform(get_req, path = tmp), error = function(e) e)
  
  # Handle 304 Not Modified case
  if (inherits(get_resp, "httr2_response") && resp_status(get_resp) == 304) {
    if (verbose) message(sprintf("[SKIP] %s not modified (304).", filename))
    unlink(tmp)
    return(invisible(local_path))
  }
  
  if (inherits(get_resp, "error")) {
    if (verbose) message(sprintf("[ERROR] GET failed for %s: %s", filename, get_resp$message))
    unlink(tmp)
    return(invisible(local_path))
  }
  
  # Success: move into place atomically
  file_move(tmp, local_path)
  
  # Capture headers
  new_etag <- resp_header(get_resp, "etag")
  new_lastmod <- resp_header(get_resp, "last-modified")
  size_bytes <- file_info(local_path)$size
  now_iso <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  
  # Upsert manifest row
  if (length(row_idx)) {
    man$etag[row_idx] <- if (!is.na(new_etag)) new_etag else remote_etag
    man$last_modified[row_idx] <- if (!is.na(new_lastmod)) new_lastmod else remote_lastmod
    man$size_bytes[row_idx] <- size_bytes
    man$downloaded_at[row_idx] <- now_iso
  } else {
    man <- rbind(man, data.frame(
      file = filename,
      etag = if (!is.na(new_etag)) new_etag else remote_etag,
      last_modified = if (!is.na(new_lastmod)) new_lastmod else remote_lastmod,
      size_bytes = size_bytes,
      downloaded_at = now_iso,
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(man, manifest_path, row.names = FALSE)
  if (verbose) message(sprintf("[OK]  %s → %s  (%,d bytes)", filename, local_path, size_bytes))
  invisible(local_path)
}

# ---- 3) Monthly updater: pass explicit files OR a regex pattern ----
# Examples:
#   bls_la_update(files = c("la.data.10.Arkansas", "la.data.10.Texas"))
#   bls_la_update(pattern = "^la\\.data\\.")
# By default, we’ll grab all la.data.* files if you set pattern accordingly.
bls_la_update <- function(files = NULL,
                          pattern = NULL,
                          verbose = TRUE) {
  if (is.null(files) && is.null(pattern)) {
    stop("Provide either `files` (vector of names) or `pattern` (regex) to select remote files.")
  }
  
  if (is.null(files)) {
    files <- bls_la_list_remote(pattern = pattern)
    if (length(files) == 0) {
      stop("No remote files matched your pattern.")
    }
  }
  
  for (f in files) {
    bls_la_download_if_new(f, verbose = verbose)
  }
  
  if (verbose) message("Done.")
}

# ---------- Examples you can run ----------

# 1) Single known file (your example)

# Big ones
  #bls_la_update(files = "la.data.0.CurrentU90-94")
  #bls_la_update(files = "la.data.0.CurrentU95-99")
  #bls_la_update(files = "la.data.0.CurrentU00-04")
  #bls_la_update(files = "la.data.0.CurrentU05-09")
  #bls_la_update(files = "la.data.0.CurrentU10-14")
  #bls_la_update(files = "la.data.0.CurrentU15-19")
  #bls_la_update(files = "la.data.0.CurrentU20-24")
  bls_la_update(files = "la.data.0.CurrentU25-29")


# Current S and U
  # bls_la_update(files = "la.data.1.CurrentS")
  # bls_la_update(files = "la.data.2.AllStatesU")
  # bls_la_update(files = "la.data.3.AllStatesS")
  # bls_la_update(files = "la.data.4.RegionDivisionU")
  # bls_la_update(files = "la.data.5.RegionDivisionS")

# Metro, Division, Micro, Combined, County, City
  # bls_la_update(files = "la.data.60.Metro")
  # bls_la_update(files = "la.data.61.Division")
  # bls_la_update(files = "la.data.62.Micro")
  # bls_la_update(files = "la.data.63.Combined")
  # bls_la_update(files = "la.data.64.County")
  # bls_la_update(files = "la.data.65.City")

# By State
  # bls_la_update(files = "la.data.51.Texas")









