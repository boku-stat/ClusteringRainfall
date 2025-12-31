get_event_characteristics <- function(start, end, x) 
{
  xx <- x %>% filter(between(date, start, end)) 
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

get_max_hourly_precip <- function(start, end, x, th)
{
  xx <- x %>% filter(between(date, start, end)) 
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

get_events <- function(fname, path_out, th_init = 0.1, nh = 4, th_sum = 1.27)
{
  x <- read_csv(fname, show_col_types = FALSE)
  x <- left_join(df, x, by = c("date" = "time"))
  x <- x %>% mutate(month = month(date))
  x <- x %>% filter(between(month, 4, 10))
  nna <- sum(is.na(x$rr))
  d <- dim(df)[1]
  if(nna > d*0.1)
  {
    print("No data")
  }
  else
  {
    x <- x %>% mutate(ev = ifelse(rr < th_init, 0, 1),
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
    ev <- x %>% group_by(event) %>% summarize(start = first(date), end = last(date)) %>% ungroup() %>% 
      filter(event != 0)
    ev <- ev %>% mutate(lag_end = lag(end, default = ev$end[1]))
    ev <- ev %>% mutate(int_period = as.integer(difftime(start, lag_end, units = "mins")))
    le <- dim(ev)[1]
    ev <- ev %>% mutate(grp_ev = 1)
    for(i in 2:le)
    {
      ev$grp_ev[i] <- ifelse(ev$int_period[i] <= nh*60, ev$grp_ev[i-1], ev$grp_ev[i-1] + 1) 
    }
    out <- ev %>% group_by(grp_ev) %>% summarize(start = min(start), end = max(end)) %>% rename(id = grp_ev)
    out <- out %>% ungroup()
    print("All events are picked")
    # Select only events that have an hourly precipitation sum larger than 1.27
    out <- out %>% mutate(th_filter = pmap_int(.l = list(start, end), .f = function(start, end)
      get_max_hourly_precip(start = start, end = end, x = x, th = th_sum)))
    out <- out %>% filter(th_filter == 1) %>% dplyr::select(!all_of(c("th_filter")))
    print("Events over threshold are filtered")
    out <- out %>% mutate(ev_char = pmap(.l = list(start, end), .f = function(start, end) 
      get_event_characteristics(start = start, end = end, x = x)))
    print("All characteristics are calculated")
    out <- unnest(out, ev_char)
    bn <- basename(fname)
    fname_out <- paste0(path_out, "events_iet_", nh, "h_", bn)
    print(fname_out)
    write_csv(out, fname_out)  
  }
}
