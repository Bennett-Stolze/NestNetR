#' Read a .lux light file
#'
#' Parses a \code{.lux} file exported from geolocators and returns a
#' two-column data frame with timestamps and log-transformed light values.
#'
#' The function looks for the first line beginning with \code{"DD/MM/YYYY"}
#' (typical header marker) and reads the tab-separated data from there.
#' Light is transformed as \eqn{\log(light + small)} and shifted to be
#' non-negative.
#'
#' @param lux_file Character scalar. Path to the \code{.lux} file.
#' @param tz Character time zone passed to \code{as.POSIXct}. Default \code{"UTC"}.
#' @param small Small positive constant added before \code{log} to avoid
#'   \code{log(0)}. Default \code{1e-4}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{Date}{POSIXct (UTC by default)}
#'   \item{Light}{non-negative, log-transformed light}
#' }
#'
#' @examples
#' \donttest{
#' # raw <- read_lux("path/to/file.lux")
#' }
#'
#' @importFrom utils read.csv 
#' @export
read_lux <- function(lux_file, tz = "UTC", small = 1e-4) {
  stopifnot(is.character(lux_file), length(lux_file) == 1)
  if (!file.exists(lux_file)) {
    stop("File not found: ", lux_file, call. = FALSE)
  }
  
  raw_lines <- base::readLines(lux_file, warn = FALSE)
  line_index <- grep("^DD/MM/YYYY", raw_lines)
  if (length(line_index) == 0) {
    stop("Could not find header marker 'DD/MM/YYYY' in file: ", lux_file, call. = FALSE)
  }
  
  dat <- utils::read.csv(
    lux_file, header = FALSE, skip = line_index,
    sep = "\t",
    col.names = c("Date", "Light"),
    colClasses = c("character", "numeric"),
    check.names = FALSE
  )
  
  # Parse time and transform light
  dat$Date  <- as.POSIXct(strptime(dat$Date, "%d/%m/%Y %H:%M:%S", tz = tz), tz = tz)
  log_light <- log(dat$Light + small)
  dat$Light <- log_light - min(log_light, na.rm = TRUE)  # shift to >= 0
  
  dat
}


#' Read a .deg file (temperature/cond data)
#'
#' Parses a \code{.deg} file exported from geolocators and returns a tidy
#' day-level data frame with minimum and maximum temperature.
#'
#' The function looks for the first line beginning with \code{"DD/MM/YYYY"}
#' (typical header marker) and reads the tab-separated data from there,
#' dropping the first three meta rows commonly present in these exports.
#'
#' @param deg_file Character scalar. Path to the \code{.deg} file.
#' @param tz Character time zone passed to \code{as.POSIXct}. Default \code{"UTC"}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{Date}{Date (day resolution)}
#'   \item{Tmin}{numeric}
#'   \item{Tmax}{numeric}
#' }
#'
#' @examples
#' \donttest{
#' # deg <- read_deg("path/to/file.deg")
#' }
#'
#' @importFrom utils read.csv 
#' @export
read_deg <- function(deg_file, tz = "UTC") {
  stopifnot(is.character(deg_file), length(deg_file) == 1)
  if (!file.exists(deg_file)) {
    stop("File not found: ", deg_file, call. = FALSE)
  }
  
  raw_lines  <- base::readLines(deg_file, warn = FALSE)
  line_index <- grep("^DD/MM/YYYY", raw_lines)
  if (length(line_index) == 0) {
    stop("Could not find header marker 'DD/MM/YYYY' in file: ", deg_file, call. = FALSE)
  }
  
  d <- utils::read.csv(
    deg_file, header = FALSE, sep = "\t", skip = line_index,
    stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE,
    col.names = c("Time", "Tmin", "Tmax", "Wets", "Conductivity")
  )
  
  # Some exports have 3 meta rows after the header marker; drop if present
  if (nrow(d) >= 3 && grepl("/", d$Time[1])) {
    d <- d[-seq_len(3), , drop = FALSE]
  }
  
  d$Date         <- as.POSIXct(strptime(d$Time, "%d/%m/%Y %H:%M:%S", tz = tz), tz = tz)
  d$Tmin         <- suppressWarnings(as.numeric(d$Tmin))
  d$Tmax         <- suppressWarnings(as.numeric(d$Tmax))
  keep           <- !is.na(d$Date) & !is.na(d$Tmin) & !is.na(d$Tmax)
  out            <- d[keep, c("Date", "Tmin", "Tmax")]
  rownames(out)  <- NULL
  out
}

#' Determine breeding period (start & end dates) from data coverage
#'
#' Uses daily coverage of light (.lux), temperature (.deg) and twilight (.twl)
#' to infer the breeding period as the longest continuous run of days with
#' light/temperature present but no twilight events.
#'
#' In `auto = TRUE` mode, the function reads the twilight file from
#' `file.path(wd, "RawData", Species, paste0(ID, "_twl.csv"))`.
#'
#' @param raw_light data.frame with POSIXct column `Date` (from \code{read_lux()}).
#' @param raw_deg   data.frame with Date or POSIXct column `Date` (from \code{read_deg()}).
#' @param ID        Character scalar; individual ID used to locate the twilight file.
#' @param Species   Character scalar; species subfolder used to locate the twilight file.
#' @param wd        Character scalar; project working directory root.
#' @param auto      Logical; if \code{TRUE} (default) read the twilight file automatically
#'                  and infer the breeding period. (Manual mode to be implemented.)
#'
#' @return Named \code{Date} vector: \code{c(start = <Date>, end = <Date>)}.
#'   If no suitable gap is found, returns \code{c(start = NA, end = NA)} with a warning.
#' @export
set_breeding_period <- function(raw_light, raw_deg, ID, Species, wd, auto = TRUE) {
  
  if (!("Date" %in% names(raw_light))) stop("raw_light must have a 'Date' column.")
  if (!("Date" %in% names(raw_deg)))   stop("raw_deg must have a 'Date' column.")
  
  if (isTRUE(auto)) {
    # ---- read twilight file ---------------------------------------------------
    twl_path <- file.path(wd, "RawData", Species, paste0(ID, "_twl.csv"))
    if (!file.exists(twl_path)) {
      stop("Twilight file not found: ", twl_path)
    }
    
    twl <- utils::read.csv(twl_path, stringsAsFactors = FALSE)
    
    # basic columns check
    if (!all(c("Twilight", "Rise", "Deleted") %in% names(twl))) {
      stop("Twilight file must contain columns: 'Twilight', 'Rise', 'Deleted'.")
    }
    
    # parse & clean
    twl$Twilight      <- as.POSIXct(twl$Twilight, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
    twl$Twilight_date <- as.Date(twl$Twilight, tz = "UTC")
    twl$Twilight_time <- format(twl$Twilight, "%H:%M:%S")
    twl                <- twl[!twl$Deleted, , drop = FALSE]
    twl                <- twl[order(twl$Twilight), , drop = FALSE]
    
    # ensure first event is sunset
    if (nrow(twl) > 0 && isTRUE(twl$Rise[1])) {
      twl <- twl[-1, , drop = FALSE]
      if (nrow(twl) == 0L) {
        warning("No twilight events left after dropping initial sunrise.")
        return(c(start = as.Date(NA), end = as.Date(NA)))
      }
    }
    
    # ---- infer breeding period from coverage ---------------------------------
    light_day <- unique(as.Date(raw_light$Date))
    deg_day   <- unique(as.Date(raw_deg$Date))
    twl_day   <- unique(as.Date(twl$Twilight_date))
    
    if (length(light_day) == 0L && length(deg_day) == 0L) {
      warning("No days available in .lux or .deg files.")
      return(c(start = as.Date(NA), end = as.Date(NA)))
    }
    
    rng_min  <- min(c(light_day, deg_day), na.rm = TRUE)
    rng_max  <- max(c(light_day, deg_day), na.rm = TRUE)
    all_days <- seq(rng_min, rng_max, by = "day")
    
    light_flag <- all_days %in% light_day
    deg_flag   <- all_days %in% deg_day
    twl_flag   <- all_days %in% twl_day
    
    # days with data but no TWL
    twl_missing <- (light_flag | deg_flag) & !twl_flag
    
    if (!any(twl_missing)) {
      warning("No gap without twilight events found within the data range.")
      return(c(start = as.Date(NA), end = as.Date(NA)))
    }
    
    r <- rle(twl_missing)
    lens_true <- ifelse(r$values, r$lengths, 0L)
    i <- which.max(lens_true)
    end_row   <- sum(r$lengths[seq_len(i)])
    start_row <- end_row - r$lengths[i] + 1L
    
    tm.breeding <- c(start = all_days[start_row], end = all_days[end_row])
  
  } else {
    # Manual path (to be implemented later)
    warning("Manual mode (auto = FALSE) not yet implemented; returning NA dates.")
    return(c(start = as.Date(NA), end = as.Date(NA)))
    
    
    tm.breeding <- ... # follows later
  }
  
  # ---- Assemble raw data for the inferred breeding period --------------------
  ## 1) Prepare 4-hourly temperature data with POSIXct Time
  deg_sub <- raw_deg[
    (if ("Time" %in% names(raw_deg)) raw_deg$Time else raw_deg$Date) >= start_buf_ct &
      (if ("Time" %in% names(raw_deg)) raw_deg$Time else raw_deg$Date) <  end_buf_ct,
    , drop = FALSE
  ]
  
  if (!("Time" %in% names(deg_sub))) {
    if (inherits(deg_sub$Date, "POSIXct")) {
      deg_sub$Time <- deg_sub$Date
    } else {
      ## If Date only is present (shouldn’t be your case now), this would create noon times:
      ## deg_sub$Time <- as.POSIXct(paste0(as.Date(deg_sub$Date), " 12:00:00"), tz="UTC")
      stop("Temperature data must have 4-hourly POSIXct timestamps in a 'Time' column.")
    }
  }
  if (!inherits(deg_sub$Time, "POSIXct")) {
    deg_sub$Time <- as.POSIXct(deg_sub$Time, tz = "UTC")
  }
  
  ## keep only the needed columns; drop rows with all-NA temps
  keep_deg <- intersect(names(deg_sub), c("Time","Tmin","Tmax"))
  deg_dat  <- deg_sub[, keep_deg, drop = FALSE]
  if (all(c("Tmin","Tmax") %in% names(deg_dat))) {
    bad <- is.na(deg_dat$Tmin) & is.na(deg_dat$Tmax)
    if (any(bad)) deg_dat <- deg_dat[!bad, , drop = FALSE]
  }
  
  ## sort both by time
  deg_dat  <- deg_dat[order(deg_dat$Time), , drop = FALSE]
  raw_light  <- raw_light[order(raw_light$Time), , drop = FALSE]
  
  ## 2) Rolling “last observation carried forward” match on time
  lux_t_num <- as.numeric(raw_light$Time)
  deg_t_num <- as.numeric(deg_dat$Time)
  
  ## for each lux time, find index of most recent deg time <= lux time
  idx <- findInterval(lux_t_num, deg_t_num)  # 0 means none earlier
  
  ## optional: enforce a maximum allowed time gap (e.g. 6 hours) to avoid carrying too far
  max_gap_sec <- 6 * 3600
  gap <- ifelse(idx > 0, lux_t_num - deg_t_num[pmax(idx, 1L)], Inf)
  idx[ gap > max_gap_sec ] <- 0L
  
  ## build output with temps from matched indices
  raw_breeding <- data.frame(
    Time  = raw_light$Time,
    Light = raw_light$Light,
    Tmin  = if ("Tmin" %in% names(deg_dat)) ifelse(idx > 0, deg_dat$Tmin[idx], NA_real_) else NA_real_,
    Tmax  = if ("Tmax" %in% names(deg_dat)) ifelse(idx > 0, deg_dat$Tmax[idx], NA_real_) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  ## optional: back-fill leading NAs (before first deg record within max_gap)
  ffill <- function(x){ last <- NA; for(i in seq_along(x)) if(!is.na(x[i])) last <- x[i] else x[i] <- last; x }
  bfill <- function(x){ nxt  <- NA; for(i in rev(seq_along(x))) if(!is.na(x[i])) nxt  <- x[i] else x[i] <- nxt;  x }
  if ("Tmin" %in% names(raw_breeding)) raw_breeding$Tmin <- bfill(ffill(raw_breeding$Tmin))
  if ("Tmax" %in% names(raw_breeding)) raw_breeding$Tmax <- bfill(ffill(raw_breeding$Tmax))
  
  ## 3) Final trim to [start, end+1 day) and drop rows without Light
  start_cut <- as.POSIXct(tm.breeding[1],      tz = "UTC")
  end_cut   <- as.POSIXct(tm.breeding[2] + 1L, tz = "UTC")  # exclusive
  raw_breeding <- raw_breeding[
    raw_breeding$Time >= start_cut &
      raw_breeding$Time <  end_cut  &
      !is.na(raw_breeding$Light),
    , drop = FALSE
  ]
  
  return(raw_breeding)
}
