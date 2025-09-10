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
#' @importFrom utils read.csv readLines
#' @export
read_lux <- function(lux_file, tz = "UTC", small = 1e-4) {
  stopifnot(is.character(lux_file), length(lux_file) == 1)
  if (!file.exists(lux_file)) {
    stop("File not found: ", lux_file, call. = FALSE)
  }
  
  raw_lines <- utils::readLines(lux_file, warn = FALSE)
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
#' @importFrom utils read.csv readLines
#' @export
read_deg <- function(deg_file, tz = "UTC") {
  stopifnot(is.character(deg_file), length(deg_file) == 1)
  if (!file.exists(deg_file)) {
    stop("File not found: ", deg_file, call. = FALSE)
  }
  
  raw_lines  <- utils::readLines(deg_file, warn = FALSE)
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
  
  d$Time         <- as.POSIXct(strptime(d$Time, "%d/%m/%Y %H:%M:%S", tz = tz), tz = tz)
  d$Tmin         <- suppressWarnings(as.numeric(d$Tmin))
  d$Tmax         <- suppressWarnings(as.numeric(d$Tmax))
  d$Date         <- as.Date(d$Time)
  keep           <- !is.na(d$Date) & !is.na(d$Tmin) & !is.na(d$Tmax)
  out            <- d[keep, c("Date", "Tmin", "Tmax")]
  rownames(out)  <- NULL
  out
}


# create segments of two days (moving)
create_segments_by_days <- function(data, date_col, segment_days = 2, overlap_days = 1) {
  # Convert date column to POSIXct if it's not already
  data <- data %>%
    mutate(!!date_col := ymd_hms(!!sym(date_col)))
  
  # Find the start and end times of the data
  min_time <- floor_date(min(data[[date_col]], na.rm = TRUE), unit = "day")
  max_time <- ceiling_date(max(data[[date_col]], na.rm = TRUE), unit = "day")
  
  # Initialize a list to store the segments
  segments <- list()
  
  # Calculate the number of days in each segment and overlap
  segment_duration <- segment_days * 24 * 60 * 60  # in seconds
  overlap_duration <- overlap_days * 24 * 60 * 60  # in seconds
  
  # Create overlapping segments
  current_start <- min_time
  segment_id <- 1
  while (current_start + segment_duration <= max_time) {
    current_end <- current_start + segment_duration
    segment <- data %>%
      filter(!!sym(date_col) >= current_start & !!sym(date_col) < current_end)
    
    segment_start_date = as.Date(current_start)  # get start date of segment
    
    segments[[segment_id]] <- segment
    segment_id <- segment_id + 1
    current_start <- current_start + overlap_duration
  }
  
  return(segments)
}

# Function to rename elements within each sublist
rename_elements_in_sublists <- function(test_data) {
  # Loop through each sublist
  for (i in seq_along(test_data)) {
    sublist <- test_data[[i]]
    
    # Loop through each element in the sublist
    for (j in seq_along(sublist)) {
      # Get the first date from the Date column of the element
      first_date <- format(sublist[[j]]$Date[1], "%Y-%m-%d")
      
      # Rename the element with its index and first date
      names(sublist)[j] <- paste0(j, "_", first_date)
    }
    
    # Update the sublist in test_data
    test_data[[i]] <- sublist
  }
  
  return(test_data)
}

# Define a recursive function to filter lists based on length (48h window)
filter_nested_list <- function(lst, target_length) {
  # If lst is a list, apply the function recursively
  if (is.list(lst)) {
    # Apply the function to each element in the list
    filtered_list <- lapply(lst, function(x) filter_nested_list(x, target_length))
    
    # Remove elements that are not of the target length
    # Keep elements that are either lists or vectors with the correct length
    filtered_list <- Filter(function(x) {
      if (is.list(x)) {
        length(x) > 0  # Ensure non-empty lists are kept
      } else if (is.vector(x)) {
        length(x) == target_length
      } else {
        FALSE
      }
    }, filtered_list)
    
    return(filtered_list)
  } else {
    # If it's not a list, check the length directly
    if (is.vector(lst) && length(lst) == target_length) {
      return(lst)
    } else {
      return(NULL)
    }
  }
}

# Rescale function 
rescale_values <- function(x, min_val, max_val) {
  scaled <- (x - min_val) / (max_val - min_val)
  # Ensure values are between 0 and 1
  scaled <- pmin(pmax(scaled, 0), 1)
  return(scaled)
}

# Function to add the list name as a column
add_list_name_and_entry_number <- function(df_list, list_name) {
  lapply(seq_along(df_list), function(i) {
    df <- df_list[[i]]
    list_segments <- strsplit(list_name, split = "_")[[1]]
    df$ID <- list_segments[1]
    df$Year <- list_segments[2]
    #df$Window <- i
    return(df)
  })
}

# Function to truncate sequences to the minimum length
truncate_sequences <- function(sequences_list, minlen) {
  truncated_sequences <- lapply(sequences_list, function(sequence) {
    sequence <- as.numeric(sequence)  # Ensure numeric type
    if (length(sequence) >= minlen) {
      sequence[1:minlen]
    }
  })
  return(do.call(rbind, truncated_sequences))  # Convert list of vectors to matrix
}
```
