#' Calculate characteristics of a single precipitation event
#'
#' Computes the event-level characteristics used as segmentation variables
#' from a 10-min precipitation time series within an event window.
#' @param start event start (datetime)
#' @param end event end (datetime)
#' @param x 10-min precipitation time series (tibble with columns
#'        \code{date} and \code{rr})
#' @return a one-row tibble with columns \code{severity} (event
#'         precipitation sum, mm), \code{magnitude} (maximum 10-min sum,
#'         mm), \code{duration} (min), \code{intensity} (severity/duration,
#'         mm min^-1), and \code{time_to_peak} (time from event start to
#'         maximum 10-min sum, expressed as a fraction of the duration)
get_event_characteristics <- function(start, end, x) 
{
  xx <- x |> filter(between(date, start, end)) 
  severity <- sum(xx$rr)
  magnitude <- max(xx$rr)
  duration <- as.integer(difftime(end, start, units = "mins")) + 10
  intensity <- severity/(duration) ### Intensity is given in mm/mins
  max_peak <- xx$date[which.max(xx$rr)]
  time_to_peak <- as.integer(difftime(max_peak, start, units = "mins")) + 10
  time_to_peak <- time_to_peak/duration
  tibble(severity = severity, magnitude = magnitude, duration = duration, intensity = intensity, 
         time_to_peak = time_to_peak)
}


#' Check whether an event exceeds an hourly precipitation threshold
#'
#' Computes the maximum hourly precipitation sum within an event (rolling
#' sum over six 10-min steps; the plain event sum for events of up to one
#' hour) and compares it against a threshold.
#' @param start event start (datetime)
#' @param end event end (datetime)
#' @param x 10-min precipitation time series (tibble with columns
#'        \code{date} and \code{rr})
#' @param th threshold for the maximum hourly precipitation sum (mm)
#' @return \code{1} if the maximum hourly sum reaches or exceeds \code{th},
#'         \code{0} otherwise
get_max_hourly_precip <- function(start, end, x, th)
{
  xx <- x |> filter(between(date, start, end)) 
  d <- as.integer(difftime(end, start, units = "mins")) + 10
  if(d <= 60)
  {
    hourly_max <- sum(xx$rr)
  }
  else
  {
    hourly_max <- max(rollsumr(xx$rr, k = 6))
  }
  return(ifelse(hourly_max >= th, 1, 0))
}

#' Separate a precipitation time series into events
#'
#' Reads a raw 10-min station time series, restricts it to the
#' April--October season, separates it into precipitation events using a
#' wet-timestep threshold and a Minimum Interevent Time (MIT) criterion,
#' filters out events below an hourly-sum threshold (see
#' \code{\link{get_max_hourly_precip}}), and attaches the event
#' characteristics (see \code{\link{get_event_characteristics}}).
#' Called for its side effect: writes
#' \code{events_iet_<nh>h_<basename(fname)>} to \code{path_out}. Skips
#' the station with a message if more than 10 \% of timesteps are missing.
#' @param fname path to the raw precipitation time series CSV (as created
#'        by \code{\link{get_csv}})
#' @param df tibble with a complete 10-min \code{date} sequence covering
#'        the study period, used to align the series
#' @param path_out output directory (character, with trailing slash)
#' @param th_init threshold above which a 10-min timestep counts as wet
#'        (mm per 10 min; default 0.1)
#' @param nh Minimum Interevent Time (h); dry gaps of up to \code{nh} hours
#'        are merged into one event (default 4)
#' @param th_sum minimum maximum-hourly precipitation sum for an event to
#'        be retained (mm; default 1.27)
get_events <- function(fname, df, path_out,
                       th_init = 0.1, nh = 4, th_sum = 1.27)
{
  x <- read_csv(fname, show_col_types = FALSE)
  x <- left_join(df, x, by = c("date" = "time"))
  x <- x |> mutate(month = month(date))
  x <- x |> filter(between(month, 4, 10))
  nna <- sum(is.na(x$rr))
  d <- dim(df)[1]
  if(nna > d*0.1)
  {
    print("No data")
  }
  else
  {
    x <- x |> mutate(ev = ifelse(rr < th_init, 0, 1),
                      lag_ev = lag(ev,1, default = 0))
    l <- dim(x)[1]
    x$event <- NA
    if(x$ev[1] == 1)
    {
      x$event[1] <- 1
    }
    else
    {
      x$event[1] <- 0
    }
    for(i in 2:l)
    {
      x$event[i] <- ifelse(x$ev[i] == 0, 0, ifelse(x$lag_ev[i] == 0, max(x$event, na.rm = TRUE) + 1, x$event[i-1]))
    }
    ev <- x |> group_by(event) |> summarize(start = first(date), end = last(date)) |> ungroup() |> 
      filter(event != 0)
    ev <- ev |> mutate(lag_end = lag(end, default = ev$end[1]))
    ev <- ev |> mutate(int_period = as.integer(difftime(start, lag_end, units = "mins")))
    le <- dim(ev)[1]
    ev <- ev |> mutate(grp_ev = 1)
    for(i in 2:le)
    {
      ev$grp_ev[i] <- ifelse(ev$int_period[i] <= nh*60, ev$grp_ev[i-1], ev$grp_ev[i-1] + 1) 
    }
    out <- ev |> group_by(grp_ev) |> summarize(start = min(start), end = max(end)) |> rename(id = grp_ev)
    out <- out |> ungroup()
    print("All events are picked")
    # Select only events that have an hourly precipitation sum larger than 1.27
    out <- out |> mutate(th_filter = pmap_int(.l = list(start, end), .f = function(start, end)
      get_max_hourly_precip(start = start, end = end, x = x, th = th_sum)))
    out <- out |> filter(th_filter == 1) |> dplyr::select(!all_of(c("th_filter")))
    print("Events over threshold are filtered")
    out <- out |> mutate(ev_char = pmap(.l = list(start, end), .f = function(start, end) 
      get_event_characteristics(start = start, end = end, x = x)))
    print("All characteristics are calculated")
    out <- unnest(out, ev_char)
    bn <- basename(fname)
    fname_out <- paste0(path_out, "events_iet_", nh, "h_", bn)
    print(fname_out)
    write_csv(out, fname_out)  
  }
}

#' Aggregate lightning data to hourly station-level flash indicators
#'
#' Filters the spatially joined lightning data to one station, aggregates
#' it to an hourly binary flash indicator, completes it against the full
#' hourly time sequence (missing hours = 0), and writes the result to
#' \code{flash_<id_sel>.csv} in \code{path_out}. Called for its side effect.
#' @param flash_data lightning observations already joined to station
#'        buffers (tibble with columns \code{id}, \code{flash}, \code{dt})
#' @param time_data tibble with a complete hourly \code{date} sequence for
#'        the study period and season
#' @param id_sel station ID to process (numeric or character, length = 1)
#' @param path_out output directory (character, with trailing slash)
save_flash <- function(flash_data, time_data,
                       id_sel, path_out)
{
  fname_out <- paste0(path_out, "flash_", id_sel, ".csv")
  xsel <- flash_data |> filter(id == id_sel)
  xsel <- xsel |> group_by(dt) |> summarize(flash = ifelse(sum(flash) > 0, 1, 0))
  out <- xsel[c("flash", "dt")] |> rename(date = dt)
  out <- left_join(time_data, out)
  out <- out |> mutate(flash = ifelse(is.na(flash), 0, flash))
  write_csv(out, fname_out)
}

#' Combine event, precipitation, and lightning data into the model input
#'
#' Matches an event-level file to its station's raw precipitation series
#' and hourly flash indicators (by station ID contained in the file
#' names), nests the raw 10-min precipitation values into each event (see
#' \code{\link{get_precip}}), and adds the binary event-level flash
#' indicator (see \code{\link{get_flash_data}}). Called for its side
#' effect: writes \code{model_input_<basename>.RDS} to \code{path_out}.
#' @param fname_event path to an event-level CSV created by
#'        \code{\link{get_events}}
#' @param path_out output directory (character, with trailing slash)
#' @param fname_flash character vector of paths to the station-level flash
#'        CSVs created by \code{\link{save_flash}}
#' @param fname_precip character vector of paths to the raw precipitation
#'        CSVs created by \code{\link{get_csv}}
combine_data <- function(fname_event, path_out,
                         fname_flash, fname_precip)
{
  bn <- basename(fname_event)
  id <- gsub(".csv", "", last(strsplit(bn, split = "_")[[1]]))
  ev <- read_csv(fname_event, show_col_types = FALSE)
  fname_raw <- grep(id, fname_precip, value = TRUE)
  x <- read_csv(fname_raw, show_col_types = FALSE)
  fname_flash <- grep(id, fname_flash, value = TRUE)
  flash <- read_csv(fname_flash, show_col_types = FALSE)
  # Add raw precipitation data to each event
  ev <- ev |> mutate(rr = pmap(.l = list(start, end), .f = function(start, end)
    get_precip(start = start, end = end, df = x)))
  ev <- ev  |>  mutate(flash = pmap_int(.l = list(start, end), .f = function(start, end) 
    get_flash_data(start = start, end = end, flash = flash)))
  fname_out <- paste0(path_out, "model_input_", bn)
  fname_out <- gsub(".csv", ".RDS", fname_out)
  saveRDS(ev, fname_out)
}