library(edfReader)
library(stringr)
library(dplyr)
library(purrr)
library(ggplot2)




# ---- Helper: read EDF header start times once (cached) ----
# Returns a named POSIXct vector: names == filenames (e.g., "chb24_01.edf")
read_edf_starts <- function(edf_dir, files, tz = "UTC") {
  unique_files <- unique(files)
  # Build full paths and read headers
  starts <- map(unique_files, function(fn) {
    path <- file.path(edf_dir, fn)
    if (!file.exists(path)) {
      warning(sprintf("EDF not found: %s", path))
      return(NA_real_)
    }
    hdr <- try(readEdfHeader(path), silent = TRUE)
    if (inherits(hdr, "try-error")) {
      warning(sprintf("Failed to read EDF header: %s", path))
      return(NA_real_)
    }
    # edfReader returns a POSIXct 'startTime' (header start date+time)
    # We force the timezone if you want all outputs in a consistent zone.
    if (!is.null(tz)) {
      return(lubridate::force_tz(hdr$startTime, tz = tz))
    } else {
      return(hdr$startTime)
    }
  })
  starts_vec <- do.call(c, starts)
  names(starts_vec) <- unique_files
  starts_vec
}

# ---- Main: parse summary + join EDF header start times ----
# txt: either a single string or character vector of lines
# edf_dir: directory containing the EDF files named in the summary
parse_seizure_text_with_edf <- function(txt, edf_dir, tz = "UTC") {
  # Accept both a character scalar or vector
  if (length(txt) == 1L) {
    con <- textConnection(txt); on.exit(close(con))
    lines <- readLines(con, warn = FALSE)
  } else {
    lines <- txt
  }
  lines <- trimws(lines)
  
  # Identify file blocks
  file_idx <- grep("^File Name:", lines)
  if (length(file_idx) == 0L) stop("No 'File Name:' entries found in summary.")
  
  out_rows <- list()
  
  for (i in seq_along(file_idx)) {
    start_line <- file_idx[i]
    end_line <- if (i < length(file_idx)) file_idx[i + 1] - 1L else length(lines)
    block <- lines[start_line:end_line]
    
    # File name
    fn_line <- block[grep("^File Name:", block)[1]]
    file_name <- sub("^File Name:\\s*", "", fn_line)
    
    # Number of seizures
    ns_pos <- grep("^Number of Seizures in File:", block)
    n_seiz <- 0L
    if (length(ns_pos)) {
      n_seiz <- suppressWarnings(as.integer(sub("^Number of Seizures in File:\\s*", "", block[ns_pos[1]])))
      if (is.na(n_seiz)) n_seiz <- 0L
    }
    if (n_seiz == 0L) next  # skip this file if no seizures
    
    # Seizure lines – support both with and without indices; allow variable spaces
    start_lines_idx <- grep("^Seizure\\s+\\d+\\s+Start Time:", block, value = TRUE)
    end_lines_idx   <- grep("^Seizure\\s+\\d+\\s+End Time:",   block, value = TRUE)
    start_lines_noi <- grep("^Seizure\\s+Start Time:", block, value = TRUE)
    end_lines_noi   <- grep("^Seizure\\s+End Time:",   block, value = TRUE)
    
    # Parsers
    get_idx <- function(x) {
      m <- str_match(x, "Seizure\\s+(\\d+)\\s+")
      idx <- suppressWarnings(as.integer(m[,2]))
      ifelse(is.na(idx), NA_integer_, idx)
    }
    get_secs <- function(x) {
      suppressWarnings(as.integer(sub(".*Time:\\s*(\\d+)\\s*seconds.*", "\\1", x)))
    }
    
    # Assemble starts/ends (indexed + unindexed)
    starts <- data.frame(idx = integer(0), secs = integer(0))
    ends   <- data.frame(idx = integer(0), secs = integer(0))
    
    if (length(start_lines_idx)) {
      starts <- rbind(starts, data.frame(
        idx  = vapply(start_lines_idx, get_idx, integer(1)),
        secs = vapply(start_lines_idx, get_secs, integer(1))
      ))
    }
    if (length(end_lines_idx)) {
      ends <- rbind(ends, data.frame(
        idx  = vapply(end_lines_idx, get_idx, integer(1)),
        secs = vapply(end_lines_idx, get_secs, integer(1))
      ))
    }
    if (length(start_lines_noi)) {
      starts <- rbind(starts, data.frame(
        idx  = seq_len(length(start_lines_noi)),
        secs = vapply(start_lines_noi, get_secs, integer(1))
      ))
    }
    if (length(end_lines_noi)) {
      ends <- rbind(ends, data.frame(
        idx  = seq_len(length(end_lines_noi)),
        secs = vapply(end_lines_noi, get_secs, integer(1))
      ))
    }
    
    # Align starts/ends by index if possible, else by order; truncate/extend using n_seiz
    if (nrow(starts) > n_seiz) starts <- starts[seq_len(n_seiz), , drop = FALSE]
    if (nrow(ends)   > n_seiz) ends   <- ends[seq_len(n_seiz), , drop = FALSE]
    
    if (nrow(starts) && nrow(ends) && all(!is.na(starts$idx)) && all(!is.na(ends$idx))) {
      se_df <- merge(starts, ends, by = "idx", all = TRUE, suffixes = c("_start", "_end"))
      se_df <- se_df[order(se_df$idx), ]
    } else {
      max_n <- max(nrow(starts), nrow(ends), n_seiz)
      se_df <- data.frame(
        idx = seq_len(max_n),
        secs_start = if (nrow(starts)) starts$secs else rep(NA_integer_, max_n),
        secs_end   = if (nrow(ends))   ends$secs   else rep(NA_integer_, max_n)
      )
    }
    
    # Emit rows (we will attach EDF header times later)
    for (k in seq_len(nrow(se_df))) {
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        file              = file_name,
        seizure_index     = se_df$idx[k],
        seizure_start_sec = se_df$secs_start[k],
        seizure_end_sec   = se_df$secs_end[k],
        stringsAsFactors  = FALSE
      )
    }
  } # next file block
  
  if (!length(out_rows)) return(data.frame())
  
  df <- bind_rows(out_rows)
  
  # Read EDF header start times once and join
  starts_vec <- read_edf_starts(edf_dir, df$file, tz = tz)  # POSIXct per file
  df <- df |>
    mutate(file_start_time = starts_vec[file],
           global_start_time = file_start_time + seizure_start_sec,
           global_end_time   = ifelse(is.na(seizure_end_sec), NA, file_start_time + seizure_end_sec))
  
  # Derive seconds (epoch and relative)
  # epoch seconds (Unix) – useful if you want raw numeric timestamps:
  df <- df |>
    mutate(
      file_start_epoch_sec   = as.numeric(file_start_time),
      global_start_epoch_sec = as.numeric(global_start_time),
      global_end_epoch_sec   = as.numeric(global_end_time)
    )
  
  # Relative seconds since the first file start (only for rows with non-NA times)
  first_start <- suppressWarnings(min(df$file_start_time, na.rm = TRUE))
  if (is.finite(as.numeric(first_start))) {
    df <- df |>
      mutate(
        global_start_sec = as.numeric(difftime(global_start_time, first_start, units = "secs")),
        global_end_sec   = as.numeric(difftime(global_end_time,   first_start, units = "secs"))
      )
  } else {
    df$global_start_sec <- NA_real_
    df$global_end_sec   <- NA_real_
  }
  
  # Order by global start (if available), else by file and index
  df <- df |>
    arrange(dplyr::desc(!is.na(global_start_sec)), global_start_sec, file, seizure_index)
  
  df
}




collapse_clusters_24 <- function(df,
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



# ---- Example usage ----
# Example: read summary text from disk and EDF headers from a directory
txt <- readLines("data/seizure-data/chb24-summary.txt")
edf_dir <- "data/seizure-data/chb24"   # directory containing chb24_*.edf

seizure_df <- parse_seizure_text_with_edf(txt, edf_dir, tz = "UTC")


ggplot(seizure_df, aes(x = global_start_sec, y = 0)) +
  geom_point(size = 3, color = "steelblue") +
  geom_rug(sides = "b") +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()



df_renewal <- collapse_clusters_24(seizure_df, min_gap_sec = 1800)


ggplot(df_renewal, aes(x = global_start_sec, y = 0)) +
  geom_point(size = 3, color = "steelblue") +
  geom_rug(sides = "b") +
  labs(
    title = "Seizure Onsets (Global Time)",
    x = "Global Time (seconds since first file)",
    y = NULL
  ) +
  theme_minimal()

head(seizure_df)
head(df_renewal)

interarrival <- diff(df_renewal$global_start_sec)   # numeric vector of interarrival seconds
hist(interarrival)

dat <- interarrival

library(marp)
set.seed(2026)


# set parameters
m <- 80 # number of iterations for MLE optimization
t <- seq(
  quantile(dat, 0.2),
  quantile(dat, 0.8),
  length.out = 4
) # time intervals
B <- 200 # number of bootstraps
BB <- 30 # number of double-bootstrapps
alpha <- 0.1 # confidence level
y <- mean(dat) 
# model_gen <- 2 # specifying the data generating model (if known)

# step one: fitting differnt renewal models
res1 <- marp::poisson_rp(dat,t,y)
res2 <- marp::gamma_rp(dat,t,m,y)
res3 <- marp::loglogis_rp(dat,t,m,y)
res4 <- marp::weibull_rp(dat,t,m,y)
res5 <- marp::lognorm_rp(dat,t,y)
res6 <- marp::bpt_rp(dat,t,m,y)

# step two: model selection and obtain model-averaged estimates
res <- marp::marp(dat,t,m,y)
res

# step three: construct different confidence intervals (including model-averaged CIs)
ci <- marp::marp_confint(dat,m,t,B,BB,alpha,y, 2)
ci


