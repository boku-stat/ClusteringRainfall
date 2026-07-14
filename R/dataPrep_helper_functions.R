#### Helper functions

#' Fit a Gamma distribution via maximum likelihood (legacy)
#'
#' Legacy helper from an earlier analysis stage; not used by the current
#' pipeline, which fits Gamma distributions via L-moments (see
#' \code{.calc_gamma()}). Requires the \code{fitdistrplus} package.
#' @param x numeric vector of (positive) precipitation values
#' @return named vector of the estimated \code{shape} and \code{rate}
#'         parameters
#' @keywords internal
get_estimate <- function(x)
{
  m <- mledist(data = x, distr = "gamma")
  m$estimate
}


#' Extract the raw precipitation values of one event
#' @param start event start (datetime)
#' @param end event end (datetime)
#' @param df 10-min precipitation time series (tibble with columns
#'        \code{time} and \code{rr})
#' @return numeric vector of 10-min precipitation sums (mm) within the
#'         event window
get_precip <- function(start, end, df)
{
  rr <- df %>% filter(between(time, start, end)) %>% pull(rr)
  rr
}

#' Derive the binary event-level lightning indicator
#'
#' Checks whether any lightning was registered between one hour before the
#' event start and the event end (hourly flash indicators, see
#' \code{\link{save_flash}}).
#' @param start event start (datetime)
#' @param end event end (datetime)
#' @param flash hourly station-level flash indicators (tibble with columns
#'        \code{date} and \code{flash})
#' @return \code{1} if lightning occurred during (or up to one hour
#'         before) the event, \code{0} otherwise
get_flash_data <- function(start, end, flash)
{
  ff <- flash %>% filter(between(date, start - hours(1), end)) %>% pull(flash)
  res <- ifelse(sum(ff) > 0, 1, 0)
  res
}