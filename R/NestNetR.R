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
#' @importFrom TwGeos lightImage
#' @importFrom dplyr filter mutate arrange
#' @importFrom lubridate hour
#' @export
set_breeding_period <- function(raw_light, raw_deg, ID, Species, wd, auto = TRUE, gr.Device = 'default') {
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
    
    # ---- Infer breeding window automatically ---------------------------------------------------
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
    return(tm.breeding)
    
  } else {
    # ---- Infer breeding window manually ---------------------------------------------------
    offset <- 0
    lmax <- 10
    zlim <- c(0,lmax)
    extend <- 0
    dark.min <- 0
    twilights <- NULL
    stage <- 1
    point.cex <- 0.8
    width <- 12
    height <- 4
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
      lightImage(raw_light, offset = offset, zlim = zlim)
      TwGeos:::selectionRectangle(Date1,Date2,col="red")
      title(main = "Select Breeding Period",
            sub = "left click: Start Point, right click: End Point, a: save period, q = quit"
            )
    }
    
    ## Draw the zoomed plot (winB)
    winBDraw <- function() {
      setDevice(winB)
      ## Display zoomed image
      lightImage(raw_light_cut, offset = offset, zlim = zlim)
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
  
#' Assemble and segment biologger data for the breeding period
#'
#' Cleans light/temperature data, subsets to breeding window, creates overlapping segments,
#' rescales values based on pre-defined quantiles, and flattens them into training-ready lists.
#'
#' @param raw_light data.frame with POSIXct-like column `Date` (from \code{read_lux()}).
#' @param raw_deg   data.frame with POSIXct-like column `Date` (from \code{read_deg()}).
#' @param tm.breeding Named vector of class \code{Date} (length 2; names 'start','end').
#' @param tz Character scalar, time zone (default "UTC").
#' @param segment_days integer, segment length in days (default 2).
#' @param overlap_days integer, overlap in days (default 1).
#' @param Light_quantiles numeric length 2, rescaling bounds.
#' @param Tmin_quantiles numeric length 2, rescaling bounds.
#' @param Tmax_quantiles numeric length 2, rescaling bounds.
#' @return A list with raw data, segmented windows, and flattened training data.
#' @importFrom magrittr %>%
#' @importFrom rlang abort warn inform
#' @importFrom dplyr filter mutate select arrange full_join any_of
#' @importFrom tidyr fill drop_na
#' @importFrom lubridate floor_date ceiling_date
#' @export
preprocessing <- function(raw_light, raw_deg, tm.breeding, tz = "UTC",
                          segment_days = 2, overlap_days = 1,
                          Light_quantiles = c(0.95, 11.1),
                          Tmin_quantiles  = c(-3.1, 22.0),
                          Tmax_quantiles  = c(2.9, 40.8)) {
  # ---- checks ----------------------------------------------------------------
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
  if (win_len < 15L)  rlang::warn(paste0("Breeding window is only ", win_len, " day(s)."))
  if (win_len > 130L) rlang::warn(paste0("Breeding window is very long (", win_len,
                                         " days) — check inputs."))
  
  # ---- normalise bounds to POSIXct in one tz ----------------------------------
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
  
  # ---- subset data -------------------------------------------------------------
  raw_light_bre <- raw_light %>%
    mutate(Date = norm_time(Date)) %>%
    filter(Date >= start_buf, Date <= end_buf)
  
  raw_deg_bre <- raw_deg %>%
    mutate(Date = norm_time(Date)) %>%
    filter(Date >= start_buf, Date <= end_buf) %>%
    select(-any_of(c("Wets", "Conductivity")))
  
  if (nrow(raw_light_bre) == 0L && nrow(raw_deg_bre) == 0L) {
    rlang::abort("No rows fall inside the (buffered) breeding window. Check dates and time zone.")
  }
  
  # ---- join & clean ------------------------------------------------------------
  raw_breeding <- raw_light_bre %>%
    full_join(raw_deg_bre, by = "Date") %>%
    arrange(Date)
  
  # duplicates check
  dup_n <- sum(duplicated(raw_breeding$Date))
  if (dup_n > 0L) {
    rlang::warn(paste0("Found ", dup_n, " duplicated timestamps after join; keeping all. ",
                       "Consider aggregating first if this is unexpected."))
  }
  
  raw_breeding <- raw_breeding %>%
    tidyr::fill(Tmin, Tmax, .direction = "up") %>%
    select(-any_of("Light_raw")) %>%
    filter(Date >= start_cut, Date < end_cut) %>%
    tidyr::drop_na()
  
  if (nrow(raw_breeding) == 0L) {
    rlang::abort("All rows dropped during cleaning. Please inspect inputs.")
  }
  
  # ---- segmentation ---------------------------------------------------------
  min_time <- lubridate::floor_date(min(raw_breeding$Date), unit = "day")
  max_time <- lubridate::ceiling_date(max(raw_breeding$Date), unit = "day")
  seg_dur  <- segment_days * 86400
  ovl_dur  <- overlap_days * 86400
  
  segments <- list(); current_start <- min_time; id <- 1
  while (current_start + seg_dur <= max_time) {
    current_end <- current_start + seg_dur
    segment <- raw_breeding %>%
      filter(Date >= current_start, Date < current_end)
    segments[[id]] <- list(
      name = paste0(id, "_", as.Date(current_start)),
      data = segment
    )
    id <- id + 1; current_start <- current_start + ovl_dur
  }
  
  # ---- rescaling ------------------------------------------------------------
  rescale <- function(x, min_val, max_val) pmin(pmax((x - min_val) / (max_val - min_val), 0), 1)
  
  segments <- lapply(segments, function(item) {
    df <- item$data %>%
      mutate(
        Light = if ("Light" %in% names(.)) rescale(Light, Light_quantiles[1], Light_quantiles[2]) else Light,
        Tmin  = if ("Tmin"  %in% names(.)) rescale(Tmin,  Tmin_quantiles[1],  Tmin_quantiles[2])  else Tmin,
        Tmax  = if ("Tmax"  %in% names(.)) rescale(Tmax,  Tmax_quantiles[1],  Tmax_quantiles[2])  else Tmax
      ) %>%
      select(Light, Tmin, Tmax)
    list(name = item$name, data = df)
  })
  
  # ---- flatten segments into vectors ----------------------------------------
  flat <- lapply(seq_along(segments), function(i) {
    item <- segments[[i]]
    parts <- strsplit(item$name, "_")[[1]]
    list(
      Light  = as.vector(item$data$Light),
      Tmin   = as.vector(item$data$Tmin),
      Tmax   = as.vector(item$data$Tmax),
      Date   = as.Date(parts[2]),
      Window = as.numeric(parts[1])
    )
  })
  
  # ---- return ---------------------------------------------------------------
  return(flat)
}

