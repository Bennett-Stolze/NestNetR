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
#' @param raw_light data.frame with POSIXct column `Date` (from \code{read_lux()}).
#' @param raw_deg   data.frame with POSIXct column `Date` (from \code{read_deg()}).
#' @param ID        Character scalar; individual ID used to locate the twilight file.
#' @param Species   Character scalar; species subfolder used to locate the twilight file.
#' @param wd        Character scalar; project working directory root.
#' @param auto      Logical; if \code{TRUE} (default), read the twilight file automatically
#'                  and infer the breeding period.
#'
#' @return A data.frame of raw biologger data restricted to the inferred breeding period,
#'   with light, Tmin, and Tmax columns.
#' @importFrom magrittr %>%
#' @importFrom dplyr filter mutate select arrange full_join left_join
#' @importFrom tidyr fill drop_na
#' @importFrom tibble tibble
#' @export
set_breeding_period <- function(raw_light, raw_deg, ID, Species, wd, auto = TRUE) {
  stopifnot("Date" %in% names(raw_light), "Date" %in% names(raw_deg))
  
  # ---- Twilight events ---------------------------------------------------------
  if (isTRUE(auto)) {
    twl_path <- file.path(wd, "RawData", Species, paste0(ID, "_twl.csv"))
    if (!file.exists(twl_path)) stop("Twilight file not found: ", twl_path)
    
    twl <- read.csv(twl_path, stringsAsFactors = FALSE) %>%
      filter(!Deleted) %>%
      mutate(
        Twilight = as.POSIXct(Twilight, tz = "UTC"),
        Twilight_date = as.Date(Twilight)
      ) %>%
      arrange(Twilight)
    
    # drop first sunrise if needed
    if (nrow(twl) > 0 && twl$Rise[1]) twl <- twl[-1, ]
  }
  
  # ---- Infer breeding window ---------------------------------------------------
  light_day <- unique(as.Date(raw_light$Date))
  deg_day   <- unique(as.Date(raw_deg$Date))
  twl_day   <- if (isTRUE(auto)) unique(as.Date(twl$Twilight_date)) else NULL
  
  all_days <- seq(min(c(light_day, deg_day), na.rm = TRUE),
                  max(c(light_day, deg_day), na.rm = TRUE),
                  by = "day")
  
  twl_missing <- (all_days %in% light_day | all_days %in% deg_day) & !(all_days %in% twl_day)
  if (!any(twl_missing)) return(c(start = NA, end = NA))
  
  r <- rle(twl_missing)
  i <- which.max(ifelse(r$values, r$lengths, 0))
  end_row   <- sum(r$lengths[seq_len(i)])
  start_row <- end_row - r$lengths[i] + 1
  tm.breeding <- c(start = all_days[start_row], end = all_days[end_row])
  
  # ---- Subset raw data ---------------------------------------------------------
  raw_light_bre <- raw_light %>%
    filter(Date >= (tm.breeding[1] - lubridate::days(1)),
           Date <= (tm.breeding[2] + lubridate::days(1))) %>%
    mutate(Date = as.POSIXct(Date, tz = "GMT"))
  
  raw_deg_bre <- raw_deg %>%
    filter(Date >= (tm.breeding[1] - lubridate::days(1)),
           Date <= (tm.breeding[2] + lubridate::days(1))) %>%
    select(-any_of(c("Wets", "Conductivity")))
  
  # ---- Join and clean ----------------------------------------------------------
  raw_breeding <- raw_light_bre %>%
    full_join(raw_deg_bre, by = "Date") %>%
    arrange(Date) %>%
    fill(Tmin, Tmax, .direction = "up") %>%
    select(-any_of("Light_raw")) %>%
    filter(Date >= tm.breeding[1],
           Date <= (tm.breeding[2] + lubridate::days(1))) %>%
    drop_na()
  
  return(raw_breeding)
}


