#' Read a light-level geolocator file
#'
#' **Reads and parses raw light intensity data** from a \code{.lux} file exported
#' from Migrate Technology \emph{Intigeo} geolocators and returns a tidy
#' two-column \code{data.frame} containing timestamps and **log-transformed,
#' normalised light values**.
#'
#' The function searches for the first line starting with \code{"DD/MM/YYYY"}
#' (typical Intigeo header marker) and reads tab-separated data from that point.
#' Light values are transformed as \eqn{\log(light +  1e-4)} to enhance low-light
#' variation and shifted to be non-negative. This transformation ensures
#' comparability with light data used for twilight estimation and deep-learning
#' input.
#'
#' @param light_file Character scalar; path to the exported light data file.
#' @param tz Character; time zone passed to \code{as.POSIXct}. Default \code{"UTC"}.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{Date}{POSIXct timestamps (UTC by default).}
#'   \item{Light}{Non-negative, log-transformed light intensity values.}
#' }
#'
#' @details
#' The resulting dataset serves as the input for further processing steps such as
#' twilight extraction (\code{TwGeos::preprocessLight}) or breeding period detection
#' within the \pkg{NestNetR} workflow.
#'
#' @examples
#' \donttest{
#' light_df <- read_light("Rawdata/D142.lux")
#' head(light_df)
#' }
#'
#' @importFrom utils read.csv
#' @export
read_light <- function(light_file, tz = "UTC") {
  stopifnot(is.character(light_file), length(light_file) == 1)
  if (!file.exists(light_file)) {
    stop("File not found: ", light_file, call. = FALSE)
  }
  
  raw_lines <- base::readLines(light_file, warn = FALSE)
  line_index <- grep("^DD/MM/YYYY", raw_lines)
  if (length(line_index) == 0) {
    stop("Could not find header marker 'DD/MM/YYYY' in file: ", light_file, call. = FALSE)
  }
  
  dat <- utils::read.csv(
    light_file, header = FALSE, skip = line_index,
    sep = "\t",
    col.names = c("Date", "Light"),
    colClasses = c("character", "numeric"),
    check.names = FALSE
  )
  
  # Parse time and transform light
  dat$Date  <- as.POSIXct(strptime(dat$Date, "%d/%m/%Y %H:%M:%S", tz = tz), tz = tz)
  log_light <- log(dat$Light + 1e-4)
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
#' @param raw_light data.frame with POSIXct column `Date` (from \code{read_light()}).
#' @param raw_deg   data.frame with POSIXct column `Date` (from \code{read_deg()}).
#' @param ID        Character scalar; individual ID used to locate the twilight file (\code{<ID>_twl.csv}).
#' @param auto      Logical; if \code{TRUE} (default), read the twilight file automatically
#'                  and infer the breeding period.
#'
#' @return A named vector containing:
#' \describe{
#'   \item{start}{Date of breeding period start (arrival at breeding site).}
#'   \item{end}{Date of breeding period end (departure from breeding site).}
#' }
#'   
#' 
#' @details
#' **Automatic detection** relies on the absence of twilight events in otherwise
#' continuous light and temperature coverage, which typically corresponds to the
#' 24-hour daylight of the Arctic breeding season. The function identifies the
#' longest such segment within the recording period and returns its start and end.
#'
#' **Manual selection** opens two interactive plots (full-year and zoomed views)
#' of the light data using \code{TwGeos::lightImage()}, allowing the user to
#' visually define the breeding window with mouse clicks:
#' \itemize{
#'   \item Left click = start of breeding period.
#'   \item Right click = end of breeding period.
#'   \item Press \code{a} to save selection, \code{q} to quit.
#' }
#' This manual mode is particularly useful when no twilight file is available or
#' when automatic detection yields ambiguous results.
#' 
#' @examples
#' \donttest{
#' # Automatic detection using twilight file
#' tm.breeding <- set_breeding_period(raw_light, raw_deg, ID = "D142", auto = TRUE)
#'
#' # Manual visual selection (e.g., when no twilight file available)
#' tm.breeding <- set_breeding_period(raw_light, raw_deg, ID = "D142", auto = FALSE)
#' }
#' 
#' @importFrom TwGeos lightImage
#' @importFrom dplyr filter mutate arrange
#' @importFrom lubridate hour month
#' @export
set_breeding_period <- function(raw_light, raw_deg, ID, auto = TRUE, gr.Device = 'default', 
                                lmax = 10, wight = 10, height = 5) {
  stopifnot("Date" %in% names(raw_light), "Date" %in% names(raw_deg))
  
  # ---- Infer breeding window automatically ---------------------------------------------------
  if (isTRUE(auto)) {
    twl_path <- file.path(dir.raw, paste0(ID, "_twl.csv"))
    if (!file.exists(twl_path)) stop("Twilight file not found: ", twl_path)
    
    twl <- read.csv(twl_path, stringsAsFactors = FALSE) |>
      filter(!Deleted) |>
      mutate(
        Twilight = as.POSIXct(Twilight, tz = "UTC"),
        Twilight_date = as.Date(Twilight)
      ) |>
      arrange(Twilight)
    
    # drop first sunrise if needed
    if (nrow(twl) > 0 && twl$Rise[1]) twl <- twl[-1, ]
    
    light_day <- unique(as.Date(raw_light$Date))
    deg_day   <- unique(as.Date(raw_deg$Date))
    twl_day   <- unique(as.Date(twl$Twilight_date))
    
    #---- Define temporal range ------------------------------------------------
    twl_start <- min(twl_day, na.rm = TRUE)
    twl_end   <- max(twl_day, na.rm = TRUE)

    # # If the last twilight date falls within the Arctic breeding season (Apr–Sep),
    # # extend trimming to the end of that calendar year to capture full breeding data
    # if (lubridate::month(twl_end) %in% 5:9) {
    #   end_of_year <- as.Date(paste0(lubridate::year(twl_end), "-12-31"))
    #   light_day <- light_day[light_day <= end_of_year]
    #   deg_day   <- deg_day[deg_day <= end_of_year]
    # } else {
      # Trim light and temperature data to the end of twilight period
      #light_day <- light_day[light_day <= twl_end]
      #deg_day   <- deg_day[deg_day <= twl_end]
    # }

    # Define complete timeline within chosen range
    all_days <- seq(twl_start, max(c(light_day, deg_day), na.rm = TRUE), by = "day")
    twl_missing <- (all_days %in% light_day | all_days %in% deg_day) & !(all_days %in% twl_day)
    if (!any(twl_missing)) return(c(start = NA, end = NA))
    
    r <- rle(twl_missing)
    i <- which(r$values)[1]
    end_row   <- sum(r$lengths[seq_len(i)])
    start_row <- end_row - r$lengths[i] + 1
    tm.breeding <- c(start = all_days[start_row], end = all_days[end_row])
    return(tm.breeding)
    
  } else {
    # ---- Infer breeding window manually ---------------------------------------------------
    zlim <- c(0,lmax)
    extend <- 0
    dark.min <- 0
    twilights <- NULL
    stage <- 1
    point.cex <- 0.8
    zoom <- 6
    raw_light_cut <- raw_light
    
    # Given a vector of POSIXct dates, extracts the time of day component of the date and returns it as decimal hours
    as.hour <- function(tm) {(as.numeric(tm)-as.numeric(as.POSIXct(as.Date(tm))))/3600}
    
    ## Round down/up to nearest offset
    floorDate <- function(date) date - ((as.hour(date)-offset)%%24)*60*60
    ceilingDate <- function(date) date + ((offset-as.hour(date))%%24)*60*60
    
    # Set starting stage & settings
    stage <- 1
    minDate <- floorDate(min(raw_light$Date))
    maxDate <- ceilingDate(max(raw_light$Date))
    Date1 <- minDate
    Date2 <- maxDate
    
    ## Set cached values
    cache <- function(k) {
      tm.breeding <<- c(as.Date(Date1), as.Date(Date2))
    }
      
    ## Select device
    setDevice <- function(dev) if(dev.cur()!=dev) dev.set(dev)
    ## Focus if possible
    focus <- setDevice
    
    # Function for mouse button clicks inside image
    ndcTsimageDate <- function(x,y) {
      day <- .POSIXct(grconvertX(x,from="ndc",to="user"),"GMT")
      hour <- grconvertY(y,from="ndc",to="user")
      .POSIXct(day+(hour%%24-hour(day))*60*60,"GMT")
    }
    
    ## Draw the logger data  (winA)
    winADraw <- function() {
      setDevice(winA)
      lightImage(raw_light, offset = 0, zlim = zlim)
      TwGeos:::selectionRectangle(Date1,Date2,col="red")
      title(main = "Select Breeding Period",
            sub = "left click: Start Point, right click: End Point, a: save period, q = quit"
            )
    }
    
    ## Draw the zoomed plot (winB)
    winBDraw <- function() {
      setDevice(winB)
      ## Display zoomed image
      lightImage(raw_light_cut, offset = 0, zlim = zlim)
      title(main = paste0("Starting Point: ", as.Date(Date1)," - ", "Ending Point: ", as.Date(Date2)))
    }
    
    ## Event handler for mouse clicks
    winAOnMouseDown <- function(buttons, x, y) {
      setDevice(winA)
      if (length(buttons) > 0) {
        b <- mouseButton(buttons)
        ## Set start- and endpoints and zoom window B
        if(b==1) {
          Date1 <<- ndcTsimageDate(x,y)
          if(Date1 > Date2) Date2 <<- maxDate
          raw_light_cut <<- raw_light[raw_light$Date >= Date1 & raw_light$Date <= Date2,]
        }
        if(b==2) {
          Date2 <<- ndcTsimageDate(x,y)
          if(Date2 < Date1) Date1 <<- minDate
          raw_light_cut <<- raw_light[raw_light$Date >= as.Date(Date1) & raw_light$Date <= as.Date(Date2),]
        }
      }
      winADraw()
      winBDraw()
      NULL
    }

    ## Function to handle mouse button presses
    mouseButton <- function(buttons) {
      n <- length(buttons)
      if(n >= 2 && buttons[1]==0 && buttons[2]==1) return(2)
      if(n >= 2 && buttons[1]==0 && buttons[2]==2) return(0)
      if(n >= 1 && buttons[1]==2) return(2)
      if(n >= 1 && buttons[1]==0) return(1)
      0
    }
    
    ## Keyboard input handler
    onKeybd <- function(key) {
      if (key == "q") return(-1)
      if (key == "a") { cache(1); return(-1) }
      winADraw()
      NULL
    }
    
    ## Set up the stages
    #if(stage == 2) cache(1)
    
    # Create the two windows
    ## winA - entire period
    if(gr.Device == 'default') dev.new(width = width, height = height, noRStudioGD = FALSE)
    if(gr.Device == "x11") x11(width = width, height = height)
    winA <- dev.cur()
    winADraw()
    setGraphicsEventHandlers(
      which = winA,
      prompt = "Select Breeding Period",
      onKeybd = onKeybd,
      onMouseDown = winAOnMouseDown)
    
    ## Create second window (winB)
    if(gr.Device=='default') dev.new(width=width,height=height,noRStudioGD=TRUE)
    if(gr.Device=="x11") x11(width=width,height=height)
    winB <- dev.cur()
    winBDraw()
    focus(winA)
    
    ## Monitor for events
    tryCatch({
      getGraphicsEvent()
      dev.off(winA)
      dev.off(winB)
    }, finally={
      dev.off(winA)
      dev.off(winB)})
    
    # bind both dates to tm.breeding
    return(tm.breeding)
  }
  return(tm.breeding)
}
  
#' Prepare breeding-period biologger data for model input
#'
#' **Assembles, cleans, and segments biologger data for the Arctic breeding period**
#' into standardised time windows suitable for deep-learning classification.  
#' Specifically, this function **subsets light and temperature data to the inferred
#' breeding window, applies rescaling based on quantile ranges, and
#' generates overlapping, fixed-length segments (windows) representing continuous
#' sequences of biologger measurements.**
#'
#' @param ID Character scalar; **individual or deployment identifier** used to tag
#'   all resulting segments.
#' @param raw_light Data frame with POSIXct column \code{Date} (from
#'   \code{read_light()}), containing **log-transformed light data**.
#' @param raw_deg Data frame with POSIXct column \code{Date} (from
#'   \code{read_deg()}), containing **minimum and maximum temperature readings**
#'   derived from the same logger.
#' @param tm.breeding Named vector of class \code{Date} (length 2; names
#'   \code{start}, \code{end}) defining the breeding period boundaries
#'   (typically returned by \code{set_breeding_period()}).
#' @param dir.breeding Character scalar; directory where processed data could be
#'   stored (currently unused but reserved for optional export as \code{.csv}).
#' @param tz Character scalar; time zone to which timestamps are converted
#'   (default \code{"UTC"}).
#' @param segment_days Integer; duration of one analysis window in days
#'   (default = **1 day**, corresponding to 288 five-minute steps).
#' @param overlap_days Integer; overlap between consecutive segments in days
#'   (default = **1 day**, i.e. fully overlapping daily windows).
#' @param Light_quantiles Numeric vector (length 2); lower and upper quantiles used
#'   for **light-value rescaling** (defaults based on the 2.5 % and 97.5 %
#'   percentiles across the Ruddy Turnstone dataset).
#' @param Tmin_quantiles Numeric vector (length 2); analogous quantiles for
#'   minimum temperature.
#' @param Tmax_quantiles Numeric vector (length 2); analogous quantiles for
#'   maximum temperature.
#'   
#' @return A **list of flattened segment objects**, each containing:
#' \describe{
#'   \item{Light, Tmin, Tmax}{Numeric vectors (length = segment_days × 288)
#'     scaled to 0–1 range.}
#'   \item{ID}{Character; deployment ID.}
#'   \item{Window}{Sequential number of the segment.}
#'   \item{Date}{Date corresponding to the segment start.}
#' 
#' @details
#' Light and temperature data are aligned, interpolated, and trimmed to the
#' breeding window before being split into overlapping windows of equal length.
#' Each variable is then rescaled using pre-defined empirical quantiles to ensure
#' comparability across individuals and seasons, matching the input expectations
#' of the CNN architecture.
#' 
#' @examples
#' \donttest{
#' tm.breeding <- set_breeding_period(raw_light, raw_deg, ID = "D142")
#' breeding_data <- preprocessing("D142", raw_light, raw_deg, tm.breeding)
#' }
#'
#' @importFrom rlang abort warn inform
#' @importFrom dplyr filter mutate select arrange full_join any_of
#' @importFrom tidyr fill drop_na
#' @importFrom lubridate floor_date ceiling_date
#' @importFrom purrr imap
#' @export
preprocessing <- function(ID, raw_light, raw_deg, tm.breeding, tz = "UTC",
                          segment_days = 1, overlap_days = 1,
                          Light_quantiles = c(0.9478656, 11.2115639),
                          Tmin_quantiles  = c(-2.8, 23.5),
                          Tmax_quantiles  = c(3.0, 40.8)) {
  # ---- checks --- 
  if (!is.data.frame(raw_light)) rlang::abort("`raw_light` must be a data.frame.")
  if (!is.data.frame(raw_deg))   rlang::abort("`raw_deg` must be a data.frame.")
  if (!("Date" %in% names(raw_light))) rlang::abort("`raw_light` needs a `Date` column.")
  if (!("Date" %in% names(raw_deg)))   rlang::abort("`raw_deg` needs a `Date` column.")
  if (!(inherits(tm.breeding, "Date") && length(tm.breeding) == 2L))
    rlang::abort("`tm.breeding` must be a Date vector of length 2 (start, end).")
  if (any(is.na(tm.breeding))) rlang::abort("`tm.breeding` contains NA values.")
  if (tm.breeding[2] < tm.breeding[1]) rlang::abort("`end` date is earlier than `start`.")
  
  
  # heuristics about window length (tweak thresholds for your species/system)
  win_len <- as.integer(diff(tm.breeding)) + 1L
  if (win_len < 40L)  rlang::warn(paste0("Breeding window for ID ", ID, " is only ", win_len, " day(s)."))
  if (win_len > 80L) rlang::warn(paste0("Breeding window for ID ", ID, " is very long (", win_len,
                                         " days) — check inputs."))
  
  # ---- normalise bounds to POSIXct in one tz ---
  # (also check if inputs came with different tz attributes)
  tz_light <- attr(raw_light$Date, "tzone")
  tz_deg   <- attr(raw_deg$Date,   "tzone")
  tzn <- unique(na.omit(c(tz_light, tz_deg)))
  if (length(tzn) > 0L && any(tzn != tz)) {
    rlang::warn(paste0("Input time zones (", paste(tzn, collapse = ", "),
                       ") differ from `tz = ", tz, "`; converting."))
  }
  
  start_buf <- as.POSIXct(as.Date(tm.breeding[1]) - 1L, tz = tz)
  end_buf   <- as.POSIXct(as.Date(tm.breeding[2]) + 1L, tz = tz)
  start_cut <- as.POSIXct(as.Date(tm.breeding[1]),     tz = tz)
  end_cut   <- as.POSIXct(as.Date(tm.breeding[2]) + 1L, tz = tz)  # end exclusive
  
  norm_time <- function(x) as.POSIXct(x, tz = tz)
  
  # ---- subset data ---
  raw_light_bre <- raw_light |>
    mutate(Date = norm_time(Date)) |>
    filter(Date >= start_buf, Date <= end_buf)
  
  raw_deg_bre <- raw_deg |>
    mutate(Date = norm_time(Date)) |>
    filter(Date >= start_buf, Date <= end_buf) |>
    select(-any_of(c("Wets", "Conductivity")))
  
  if (nrow(raw_light_bre) == 0L && nrow(raw_deg_bre) == 0L) {
    rlang::abort("No rows fall inside the (buffered) breeding window. Check dates and time zone.")
  }
  
  # ---- join & clean ---
  raw_breeding <- raw_light_bre |>
    full_join(raw_deg_bre, by = "Date") |>
    arrange(Date)
  
  # duplicates check
  dup_n <- sum(duplicated(raw_breeding$Date))
  if (dup_n > 0L) {
    rlang::warn(paste0("Found ", dup_n, " duplicated timestamps after join; keeping all. ",
                       "Consider aggregating first if this is unexpected."))
  }
  
  raw_breeding <- raw_breeding |>
    tidyr::fill(Tmin, Tmax, .direction = "up") |>
    select(-any_of("Light_raw")) |>
    filter(Date >= start_cut, Date < end_cut) |>
    tidyr::drop_na()
  
  if (nrow(raw_breeding) == 0L) {
    rlang::abort("All rows dropped during cleaning. Please inspect inputs.")
  }
  
  # ---- segmentation ---
  min_time <- lubridate::floor_date(min(raw_breeding$Date), unit = "day")
  max_time <- lubridate::ceiling_date(max(raw_breeding$Date), unit = "day")
  seg_dur  <- segment_days * 86400
  ovl_dur  <- overlap_days * 86400
  
  segments <- list(); current_start <- min_time; id <- 1
  while (current_start + seg_dur <= max_time) {
    current_end <- current_start + seg_dur
    
    segment <- raw_breeding |>
      filter(Date >= current_start, Date < current_end)
    
    # Create the segment name directly
    seg_name <- paste0(id, "_", as.Date(current_start))
    
    # Assign the data.frame directly to that list name
    segments[[seg_name]] <- segment
    
    id <- id + 1
    current_start <- current_start + ovl_dur
  }
  
  # ---- filter incomplete segments ---
  expected_len <- segment_days * 24 * (60 / 5)
  segments <- segments[sapply(segments, nrow) == expected_len]
  
  if (length(segments) == 0L) {
    rlang::warn("No complete segments retained after filtering — check input resolution or window size.")
  }
  
  # ---- rescaling ---
  rescale <- function(x, min_val, max_val) pmin(pmax((x - min_val) / (max_val - min_val), 0), 1)
  
  segments <- lapply(segments, function(df) {
    df |>
      mutate(
        Light = if ("Light" %in% names(df)) rescale(Light, Light_quantiles[1], Light_quantiles[2]) else Light,
        Tmin  = if ("Tmin"  %in% names(df)) rescale(Tmin,  Tmin_quantiles[1],  Tmin_quantiles[2])  else Tmin,
        Tmax  = if ("Tmax"  %in% names(df)) rescale(Tmax,  Tmax_quantiles[1],  Tmax_quantiles[2])  else Tmax
      ) |>
      select(Light, Tmin, Tmax)
  })
  
  # ---- flatten segments into vectors ---
  flat <- imap(segments, function(seg, list_name) {
    seg <- as.data.frame(seg)
    
    dat <- seg |>
      mutate(
        Light = if ("Light" %in% names(seg)) as.numeric(Light) else NULL,
        Tmin  = if ("Tmin"  %in% names(seg)) as.numeric(Tmin)  else NULL,
        Tmax  = if ("Tmax"  %in% names(seg)) as.numeric(Tmax)  else NULL
      ) |>
      select(any_of(c("Light", "Tmin", "Tmax")))
    
    parts <- strsplit(list_name, "_")[[1]]
    
    list(
      Light  = if ("Light" %in% names(dat)) as.vector(dat$Light) else NA,
      Tmin   = if ("Tmin"  %in% names(dat)) as.vector(dat$Tmin)  else NA,
      Tmax   = if ("Tmax"  %in% names(dat)) as.vector(dat$Tmax)  else NA,
      ID     = ID,
      Window = as.numeric(parts[1]),
      Date   = as.Date(parts[2])
    )
  })
  return(flat)
}

#' Classify breeding behaviour from segmented biologger data
#'
#' Applies a pre-trained deep learning model to predict breeding behaviours
#' (incubation, brooding, random) for each segmented biologger window.
#' The function loads a pre-trained CNN model, generates predictions from
#' light- and temperature-based input sequences, and formats results into a
#' clean data frame ready for ecological analysis or visualisation.
#'
#' @param breeding_data A list of segmented and normalised input data, typically
#'   returned by \code{preprocessing()} or similar, where each list element
#'   contains \code{Light}, \code{Tmin}, and \code{Tmax} vectors, plus metadata
#'   fields (\code{ID}, \code{Window}, \code{Date}).
#' @param model Character; either the keyword \code{"base"} (default, loading
#'   \code{Model/base_model.keras} from the working directory) or a full path to a
#'   custom \code{.keras} file.
#'
#' @details
#' Internally, the function truncates all input channels to a common window
#' length, combines them into a 3-channel tensor, loads the pre-trained model
#' from \code{<working directory>/Model/base_model.keras}, and performs class
#' predictions. Probabilities are returned for all defined classes, together
#' with the most likely predicted class per segment.
#'
#' @return
#' A data frame with one row per segment and columns:
#' \describe{
#'   \item{ID}{Individual or deployment ID.}
#'   \item{Window}{Segment index.}
#'   \item{Date}{Start date of the segment.}
#'   \item{Year}{Calendar year extracted from \code{Date}.}
#'   \item{brooding, incubation, random}{Predicted class probabilities.}
#'   \item{Class}{Predicted behavioural class label.}
#' }
#'
#' @examples
#' \dontrun{
#' predictions <- classify_breeding_behaviour(breeding_data)
#' head(predictions)
#' }
#'
#' @importFrom keras3 load_model
#' @importFrom abind abind
#' @importFrom dplyr mutate relocate
#' @importFrom lubridate year
#' @export
classify_breeding_behaviour <- function(breeding_data, model = "base") {
  # --- find minimum sequence length automatically ---
  seq_lengths <- sapply(breeding_data, function(x) length(x$Light))
  minlen <- min(seq_lengths, na.rm = TRUE)
  
  # --- preprocess for model input ---
  light_truncated <- do.call(rbind, lapply(breeding_data, function(x) {
    vals <- as.numeric(x$Light)
    if (length(vals) >= minlen) vals[1:minlen]
  }))
  
  tmin_truncated <- do.call(rbind, lapply(breeding_data, function(x) {
    vals <- as.numeric(x$Tmin)
    if (length(vals) >= minlen) vals[1:minlen]
  }))
  
  tmax_truncated <- do.call(rbind, lapply(breeding_data, function(x) {
    vals <- as.numeric(x$Tmax)
    if (length(vals) >= minlen) vals[1:minlen]
  }))
  x_breeding <- abind::abind(light_truncated, tmin_truncated, tmax_truncated, along = 3)
  
  # --- load model ---
  if (model == "base") {
    model <- file.path(wd, "Model", "base_model.keras")
  } 
  
  if (grepl("\\.keras$", model) == FALSE) {
    stop("Invalid model specified. Use 'base' or indicate correct file with .keras extension.")
  }
  model <- keras3::load_model(model)
  
  # --- Predictions ---
  predictions <- model |> predict(x_breeding)
  
  # --- process prediction output ---
  # Convert predictions to a data frame
  classes <- c("brooding", "incubation", "random")
  predictions <- as.data.frame(predictions)
  colnames(predictions) <- levels(factor(classes))
  
  # Get the predicted classes
  predicted_classes <- apply(predictions, 1, which.max)
  predicted_classes <- levels(factor(classes))[predicted_classes]
  
  # --- reassemble metadata ---
  breeding_data_df <- lapply(breeding_data, function(x) {
    x$Light <- NULL  # Remove the "Light" entry from each list
    x$Light_raw <- NULL
    x$Tmin <- NULL
    x$Tmax <- NULL
    x$Wets <- NULL
    x$Conductivity <- NULL
    x$Class <- NULL
    x$Label <- NULL
    as.data.frame(x)  # Convert the modified list to a data frame
    return(x)  # Return the modified list
  })
  
  # --- combine everything ---
  breeding_data_df <- do.call(rbind, breeding_data_df)
  
  predictions <- cbind(breeding_data_df, predictions)
  predictions$Class <- predicted_classes
  rownames(predictions) <- NULL
  
  predictions <- predictions |>
    dplyr::mutate(Date = as.Date(unlist(Date)),
                  Year = lubridate::year(Date)) |>
    dplyr::relocate(Year, .before = Window)
  
  return(predictions)
}

#' Create list of segmented breeding period datasets for mulitple individuals
#'
#' **Automates the full preprocessing workflow across multiple geolocator files.**
#' Reads individual light and temperature files from a directory, determines each
#' bird’s breeding period, preprocesses the data, and stores the resulting
#' segmented time series per individual. Each list element corresponds to one
#' individual ID and contains a list of its breeding-period segments.
#'
#' @param dir.raw Character; directory containing raw light and temperature files
#'   (e.g. `"data/RawData/SpeciesName/"`).
#' @param auto Logical; if `TRUE` (default), breeding periods are detected
#'   automatically using \code{set_breeding_period()}. Otherwise, manually annotation is required.
#' @param segment_days Integer; length of each segmentation window in days
#'   (default = **1 day**, corresponding to 288 five-minute steps).
#'
#' @return A **named list** where each element corresponds to one individual ID
#'   and contains that individual’s segmented and normalised breeding-period data.
#'   Each sub-element is itself a list of window objects with numeric vectors
#'   (\code{Light}, \code{Tmin}, \code{Tmax}) and metadata fields (\code{ID},
#'   \code{Window}, \code{Date}).
#'
#' @examples
#' \donttest{
#' breeding_data_list <- create_breeding_data_list("data/RawData/RuddyTurnstone/", auto = TRUE)
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export
create_breeding_data_list <- function(dir.raw, auto = TRUE, segment_days = 1) {
  ID_list <- sub("\\.lux$", "", basename(list.files(dir.raw, pattern = ".lux", full.names = TRUE)))
  Species <- strsplit(dir.raw, "/")[[1]][[3]]
  
  # Create emtpy list
  breeding_data_list <- list()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(ID_list), style = 3, width = 25)
  
  for (ID in ID_list) {
    raw_light <- NestNetR::read_light(file.path(dir.raw, paste0(ID, ".lux")))
    raw_deg <- NestNetR::read_deg(file.path(dir.raw, paste0(ID, ".deg")))
    
    tm.breeding <- NestNetR::set_breeding_period(raw_light, raw_deg, ID, auto)
    
    breeding_data <- NestNetR::preprocessing(ID, raw_light, raw_deg, tm.breeding, segment_days = segment_days)
    
    breeding_data_list[[ID]] <- breeding_data
    
    setTxtProgressBar(pb, which(ID_list == ID))
  }
  
  close(pb)
  return(breeding_data_list)
}

#' Interactively label biologger segments for breeding-behaviour training data
#'
#' Opens an interactive graphical interface to manually assign behavioural classes
#' (*incubation*, *brooding*, *random*) to randomly selected segments
#' from a list of preprocessed biologger data. This function is designed to support
#' the manual creation of balanced training data for the deep-learning model.
#'
#' @param breeding_data_list A named list of segmented breeding-period data,
#'   typically created by \code{create_breeding_data_list()}, where each list element
#'   corresponds to one segment (e.g. \code{D142_1_2013-05-20}) containing light and
#'   temperature time series as numeric vectors.
#' @param segment_days Integer; the duration of one segment in days (used to display
#'   correct time bounds in the plot).
#' @param lmax Numeric; maximum light intensity used for visual scaling (default = 1).
#' @param width,height Numeric; dimensions of the plotting window in inches.
#' @param gr.Device Character; the graphics device used for interactive display.
#'   Accepts \code{"default"} (RStudioGD) or \code{"x11"} for external windows.
#'
#' @details
#' The function displays a random unlabelled segment in an interactive plot, where
#' the user can inspect light and temperature data and assign behavioural categories
#' using keyboard shortcuts:
#' \itemize{
#'   \item \code{b} – mark as *brooding*
#'   \item \code{i} – mark as *incubation*
#'   \item \code{r} – mark as *random*
#'   \item \code{s} – skip current sample
#'   \item \code{u} – undo last classification
#'   \item \code{q} – quit labelling session
#' }
#'
#' The interface also displays the current counts of each class to ensure
#' approximately balanced sample sizes. Each labelled segment is updated in the
#' \code{breeding_data_list} object by adding the fields \code{Label = "training"}
#' and \code{Class = <behaviour>}.
#'
#' @return
#' The updated \code{breeding_data_list} including manual class labels for all
#' annotated segments.
#'
#' @note
#' This function requires an interactive graphics device (it cannot be run in
#' non-interactive or headless environments).
#'
#' @examples
#' \dontrun{
#' # Start interactive training data selection
#' breeding_data_list <- create_trainingdata(
#'   breeding_data_list,
#'   segment_days = 1,
#'   gr.Device = "x11"
#' )
#' }
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom lubridate days
#' @importFrom TwGeos lightImage
#' @export
create_trainingdata <- function(breeding_data_list, segment_days, zlim=c(0,10), width=10, height=5, gr.Device='default'){
  random_count <- brooding_count <- incubation_count <- 0
  history <- list() 
  quit <- FALSE
  
  id_summary <- do.call(rbind, lapply(names(breeding_data_list), function(id) {
    segs <- breeding_data_list[[id]]
    
    # collect all Date vectors while preserving class
    all_dates <- do.call(c, lapply(segs, function(x) x$Date))
    
    data.frame(
      ID         = id,
      start_date = min(all_dates, na.rm = TRUE),
      end_date   = max(all_dates, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  breeding_data_list <- unlist(breeding_data_list, recursive = FALSE)
  names(breeding_data_list) <- gsub("\\.", "_", names(breeding_data_list))
  
  repeat{
    # Create a list of candidates from breeding_data (which are not yet assigned for "training")
    candidates <- names(breeding_data_list)[
      sapply(breeding_data_list, function(x) {
        !("Label" %in% names(x)) || is.null(x$Label) || all(is.na(x$Label))
      })
    ]
    
    if(length(candidates)==0){ message("No more candidates."); break }
    
    chosen <- sample(candidates,1) 
    parts <- strsplit(chosen,"_")[[1]]
    ID <- parts[1]
    window <- parts[2]
    Date1 <- as.POSIXct(paste0(parts[3]," 00:00:00"), tz="UTC")
    Date2 <- Date1 + segment_days*24*60*60
    start_date <- id_summary[id_summary$ID==ID,"start_date"]
    end_date <- id_summary[id_summary$ID==ID,"end_date"]
    tagdata <- subset(read_light(file.path(wd,"RawData",Species,paste0(ID,".lux"))),
                      Date >= start_date & Date <= end_date)
    
    # --- ensure Date1 and Date2 fall within available data ---
    if (nrow(tagdata) == 0) next  # skip empty tagdata
    
    data_min <- min(tagdata$Date, na.rm = TRUE)
    data_max <- max(tagdata$Date, na.rm = TRUE)
    
    # If the selected window starts before or after available data, adjust it
    if (Date1 < data_min) {
      Date1 <- data_min
      Date2 <- Date1 + lubridate::days(segment_days)
    }
    if (Date2 > data_max) {
      Date2 <- data_max
      Date1 <- Date2 - lubridate::days(segment_days)
    }
    
    # Edge case: if still out of range or invalid after adjustment, skip
    if (Date1 < data_min | Date2 > data_max | Date1 >= Date2) next
    
    behaviour_label <- NULL 
    skip <- FALSE 
    undo <- FALSE
    
    ## Set cached values
    cache <- function(k) {
      # Get label choice
      behaviour_label <<- behaviour_label
    }
    ## Select device
    setDevice <- function(dev) if(dev.cur()!=dev) dev.set(dev)
    
    ## Draw the behaviour identification window (winA)
    winADraw <- function(){
      setDevice(winA); lightImage(tagdata, zlim=zlim, offset=0)
      TwGeos:::selectionRectangle(Date1,Date2,col="red"); rect(Date1,0.1,Date2,23.9,border="red",lwd=4)
      title(main="Define behaviour", sub="(b=brooding, r=random, i=incubation, s=skip, u=undo, q=quit)")
      mtext(paste0("Brooding: ",brooding_count,", Incubation: ",incubation_count,", Random: ",random_count), cex=.8)
      mtext(paste0("ID: ",ID," | Window: ",window," | Date: ",format(Date1,"%Y-%m-%d")), side=3, line=-1.5, cex=.9, col="darkred", font=2)
    }
    onKeybd <- function(key){
      if(key=="q"){ quit <<- TRUE; return(-1) }
      if(key=="s"){ skip <<- TRUE; return(-1) }
      if(key=="u"){ undo <<- TRUE; return(-1) }
      if(key=="i"){ behaviour_label <<- "incubation"; cache(1); return(-1) }
      if(key=="r"){ behaviour_label <<- "random";     cache(1); return(-1) }
      if(key=="b"){ behaviour_label <<- "brooding";   cache(1); return(-1) }
      winADraw(); NULL
    }
    
    if(gr.Device=='default') dev.new(width=width, height=height, noRStudioGD=FALSE) else if(gr.Device=="x11") x11(width=width, height=height)
    winA <- dev.cur(); winADraw()
    tryCatch({ setGraphicsEventHandlers(which=winA, prompt="Create Training Data", onKeybd=onKeybd); getGraphicsEvent(); dev.off(winA) },
             finally = { if(dev.cur()==winA) dev.off(winA) })
    
    if(quit) break
    if(skip) next
    
    # --- Undo the last label ---
    if (undo && length(history) > 0) {
      last <- history[[length(history)]]
      seg_name <- last$seg  # full segment name, e.g. "D148_1_2013-05-20"
      
      if (seg_name %in% names(breeding_data_list)) {
        seg <- breeding_data_list[[seg_name]]
        
        # Remove labels (reset segment)
        seg$Label <- NULL
        seg$Class <- NULL
        
        breeding_data_list[[seg_name]] <- seg
        
        # Adjust class counters
        if (last$class == "incubation") incubation_count <- incubation_count - 1
        if (last$class == "brooding")   brooding_count   <- brooding_count - 1
        if (last$class == "random")     random_count     <- random_count - 1
        
        # Remove last action from history
        history <- history[-length(history)]
        
        message("Last label undone for: ", seg_name)
      } else {
        warning("Segment not found in breeding_data_list: ", seg_name)
      }
      
      next  # Continue to next iteration
    }
    
    # apply label (as attributes, not columns)
    if (exists("behaviour_label")) {
      if (behaviour_label == "incubation") {
        # Mark this sample as used in breeding_data
        breeding_data_list[[chosen]]$Label <- "training"
        # Add behaviour label to chosen sample
        breeding_data_list[[chosen]]$Class <- behaviour_label
        # set count of behaviour + 1
        incubation_count <- incubation_count + 1
      } else if (behaviour_label == "brooding") {
        # Mark this sample as used in breeding_data
        breeding_data_list[[chosen]]$Label <- "training"
        # Add behaviour label to chosen sample
        breeding_data_list[[chosen]]$Class <- behaviour_label
        # set count of behaviour + 1
        brooding_count <- brooding_count + 1
      } else if (behaviour_label == "random") {
        # Mark this sample as used in breeding_data
        breeding_data_list[[chosen]]$Label <- "training"
        # Add behaviour label to chosen sample
        breeding_data_list[[chosen]]$Class <- behaviour_label
        # set count of behaviour + 1
        random_count <- random_count + 1
      }
      # save in history
      history[[length(history)+1]] <- list(chosen=chosen, class=behaviour_label)
    }
    breeding_data_list
  }
  preclassified_data <- breeding_data_list[
    sapply(breeding_data_list, function(x) {
      "Label" %in% names(x) && !is.null(x$Label) && identical(x$Label, "training")
    })]
  preclassified_data
}
