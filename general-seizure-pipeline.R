library(stringr)
library(lubridate)
library(dplyr)
library(edfReader)


# Parse seizure summary text into a data frame using arbitrary date
parse_seizure_text <- function(txt, origin_date = "1970-01-01", tz = "UTC") {
  # Accept either a single character string or a character vector (lines)
  if (length(txt) == 1L) {
    con <- textConnection(txt)
    on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
  } else {
    lines <- txt
  }
  lines <- trimws(lines)
  
  # Helper: parse "H:M:S" or "HH:MM:SS" into seconds-of-day; supports hour 0-24 (24 -> 0 + next day)
  hms_to_sod <- function(hms) {
    # Remove any stray spaces
    hms <- gsub("\\s+", "", hms)
    parts <- strsplit(hms, ":", fixed = TRUE)[[1]]
    if (length(parts) != 3) stop(paste("Invalid time format:", hms))
    h <- as.integer(parts[1]); m <- as.integer(parts[2]); s <- as.integer(parts[3])
    if (is.na(h) || is.na(m) || is.na(s)) stop(paste("Invalid time components:", hms))
    if (h == 24L) {
      # Represent 24:MM:SS as 0:MM:SS and indicate caller should handle day rollover if needed.
      return(list(sod = m * 60 + s, hour24 = TRUE))
    } else {
      return(list(sod = h * 3600 + m * 60 + s, hour24 = FALSE))
    }
  }
  
  # Find file block starts
  file_idx <- grep("^File Name:", lines)
  if (length(file_idx) == 0L) stop("No 'File Name:' entries found.")
  
  origin <- as.POSIXct(paste0(origin_date, " 00:00:00"), tz = tz)
  day_offset <- 0L
  prev_start_sod <- NA_integer_
  first_global_file_start <- NA
  out_rows <- list()
  
  for (i in seq_along(file_idx)) {
    start_line <- file_idx[i]
    end_line <- if (i < length(file_idx)) file_idx[i + 1] - 1L else length(lines)
    block <- lines[start_line:end_line]
    
    # Extract file name
    fn_line <- block[grep("^File Name:", block)[1]]
    file_name <- sub("^File Name:\\s*", "", fn_line)
    
    # Extract file start time
    fs_pos <- grep("^File Start Time:", block)
    if (length(fs_pos) == 0L) next
    fs_line <- block[fs_pos[1]]
    fs_str <- sub("^File Start Time:\\s*", "", fs_line)
    fs_parsed <- hms_to_sod(fs_str)
    file_start_sod <- fs_parsed$sod
    # Increment day offset if time-of-day goes backwards (cross midnight)
    if (!is.na(prev_start_sod) && file_start_sod < prev_start_sod) {
      day_offset <- day_offset + 1L
    }
    prev_start_sod <- file_start_sod
    
    # Global file start POSIXct
    file_global_start <- origin + day_offset * 86400 + file_start_sod
    if (is.na(first_global_file_start)) first_global_file_start <- file_global_start
    
    # Extract number of seizures
    ns_pos <- grep("^Number of Seizures in File:", block)
    n_seiz <- 0L
    if (length(ns_pos)) {
      n_seiz <- as.integer(sub("^Number of Seizures in File:\\s*", "", block[ns_pos[1]]))
      if (is.na(n_seiz)) n_seiz <- 0L
    }
    if (n_seiz == 0L) next  # skip files without seizures
    
    # Extract seizure start/end seconds (handle variable spacing)
    start_lines <- grep("^Seizure\\s+\\d+\\s+Start Time:", block, value = TRUE)
    end_lines   <- grep("^Seizure\\s+\\d+\\s+End Time:",   block, value = TRUE)
    
    # Build a map from seizure index -> start/end seconds
    # Robust to inconsistent spaces and extra spaces before 'seconds'
    get_idx <- function(x) as.integer(sub(".*Seizure\\s+(\\d+)\\s+.*", "\\1", x))
    get_secs <- function(x) as.integer(sub(".*Time:\\s*(\\d+)\\s*seconds.*", "\\1", x))
    
    starts <- data.frame(
      idx  = vapply(start_lines, get_idx, integer(1)),
      secs = vapply(start_lines, get_secs, integer(1)),
      stringsAsFactors = FALSE
    )
    ends <- data.frame(
      idx  = vapply(end_lines, get_idx, integer(1)),
      secs = vapply(end_lines, get_secs, integer(1)),
      stringsAsFactors = FALSE
    )
    # Merge by seizure index if possible; if mismatch, record what's available
    se_df <- merge(starts, ends, by = "idx", all = TRUE, suffixes = c("_start", "_end"))
    se_df <- se_df[order(se_df$idx), ]
    
    # Create rows
    for (k in seq_len(nrow(se_df))) {
      sz_idx <- se_df$idx[k]
      sz_start_sec <- se_df$secs_start[k]
      sz_end_sec   <- se_df$secs_end[k]
      
      # Compute global start/end
      global_start_time <- file_global_start + sz_start_sec
      global_end_time   <- if (!is.na(sz_end_sec)) file_global_start + sz_end_sec else NA
      global_start_sec  <- as.numeric(difftime(global_start_time, first_global_file_start, units = "secs"))
      global_end_sec    <- if (!is.na(global_end_time)) as.numeric(difftime(global_end_time, first_global_file_start, units = "secs")) else NA
      
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        file                = file_name,
        seizure_index       = sz_idx,
        seizure_start_sec   = sz_start_sec,
        seizure_end_sec     = sz_end_sec,
        file_start_time     = file_global_start,
        # Optional: keep file time-of-day for reference
        file_start_hms      = fs_str,
        global_start_sec    = global_start_sec,
        global_end_sec      = global_end_sec,
        global_start_time   = global_start_time,
        global_end_time     = global_end_time,
        stringsAsFactors    = FALSE
      )
    }
  }
  
  if (length(out_rows) == 0L) {
    return(data.frame())
  }
  df <- do.call(rbind, out_rows)
  
  # Order by global start time to be safe
  df <- df[order(df$global_start_sec), ]
  rownames(df) <- NULL
  df
}


collapse_clusters <- function(df,
                              min_gap_sec,
                              within_file = TRUE) {
  stopifnot(all(c("file", "global_start_sec") %in% names(df)))
  if (!("global_end_sec" %in% names(df))) df$global_end_sec <- NA_real_
  
  # Order deterministically
  df <- df |>
    arrange(file, global_start_sec, seizure_index)
  
  # Choose grouping key
  grp <- if (within_file) "file" else NULL
  
  df_tag <- df |>
    group_by(across(all_of(grp))) |>
    arrange(global_start_sec, .by_group = TRUE) |>
    mutate(
      prev_end   = lag(global_end_sec),
      prev_start = lag(global_start_sec),
      # Gap is start - previous end; if previous end missing, use previous start
      gap = global_start_sec - if_else(is.na(prev_end), prev_start, prev_end),
      # New cluster when no previous event or gap >= min_gap_sec
      new_event = if_else(is.na(gap) | gap >= min_gap_sec, 1L, 0L),
      event_id  = cumsum(new_event)
    ) |>
    ungroup()
  
  # Summarise per cluster
  out <- df_tag |>
    group_by(!!!rlang::syms(c(grp, "event_id"))) |>
    summarise(
      # Representative onset = earliest in cluster
      global_start_sec = min(global_start_sec, na.rm = TRUE),
      # Cluster end: if all ends are NA, fallback to last start
      global_end_sec   = {
        all_na_end <- all(is.na(global_end_sec))
        if (all_na_end) max(global_start_sec, na.rm = TRUE)
        else max(global_end_sec, na.rm = TRUE)
      },
      n_seizures       = n(),
      .groups = "drop"
    ) |>
    arrange(!!!rlang::syms(grp), global_start_sec, event_id)
  
  out
}
