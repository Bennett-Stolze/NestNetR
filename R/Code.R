# ========================================================================== #
# Read .lux/.lig files ====
# ========================================================================== #
read_lux <- function(lux_file) {
  # Investigate the line to start reading in the data
  raw_lux <- readLines(lux_file)
  line_index <- grep("^DD/MM/YYYY", raw_lux)
  
  # Light Data loading and preprocessing
  raw_lux <- read.csv(lux_file, header = FALSE, skip = line_index, sep = "\t", 
                      col.names = c("Date", "Light"), colClasses = c("character", 
                                                                     "numeric"))
  raw_lux$Date <- as.POSIXct(strptime(raw_lux$Date, "%d/%m/%Y %H:%M:%S", 
                                      tz = "GMT"))
  names(raw_lux) <- c("Date", "Light")
  raw_lux$Light <- log(raw_lux$Light + 0.0001) + abs(min(log(raw_lux$Light + 0.0001)))
  raw_lux$Date <- as.POSIXct(raw_lux$Date, tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
    
  return(raw_lux)
}


# ========================================================================== #
# Create .deg file ====
# ========================================================================== #
read_deg <- function(deg_file) {
  # investigate the line to start reading in the data
  d <- readLines(deg_file)
  line_index <- grep("^DD/MM/YYYY", d)
  
  d <- read.csv(deg_file, header = F, sep = "\t", skip = line_index,
                stringsAsFactors = FALSE, check.names = FALSE, fill = TRUE, 
                col.names = c("Time", "Tmin", "Tmax", "Wets", "Conductivity"))[-c(1:3),]
  d$Time <- as.POSIXct(strptime(d$Time, "%d/%m/%Y %H:%M:%S", 
                                tz = "GMT"))
  d$Tmin <- as.numeric(d$Tmin)
  d$Tmax <- as.numeric(d$Tmax)
  d$Wets <- as.numeric(d$Wets)
  d$Conductivity <- as.numeric(d$Conductivity)
  d$Date <- as.Date(d$Time, format = "%d/%m/%Y")
  d <- d %>% 
    select(Date, Tmin, Tmax)
  
  if (any(is.na(d))) {
    # delete rows with NA values
    d <- d %>% drop_na()
  }
  return(d)
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
