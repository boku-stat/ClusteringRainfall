
# Source: N/nwnr/Projekte/EROSA-Stat/Tawes_stations/download_station_data.R, 01.10.25
# split into functions and script (<-- in inst/001_download_station_data.R)


################################################################################
########## Functions for downloading TAWES stations data #######################
################################################################################

#' Build a Geosphere API query URL for 10-min station data
#'
#' Constructs the download URL for the Geosphere Austria "klima-v2-10min"
#' historical station dataset, returning results in CSV format.
#' @param id station ID (numeric or character, length = 1)
#' @param valid_from start of the requested time range (ISO 8601 datetime string)
#' @param valid_to end of the requested time range (ISO 8601 datetime string)
#' @param parameter name of the climatic variable to download (character,
#'        length = 1; \code{"rr"} for 10-min precipitation sums)
#' @return the query URL (character, length = 1)
get_url <- function(id, valid_from, valid_to, parameter)
{
  url <- paste0("https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min?parameters=", parameter, "&start=", valid_from, 
                "&end=", valid_to, "&station_ids=", id, "&output_format=csv")
  url
}

#' Download a station time series and write it to disk
#'
#' Reads the CSV returned by a Geosphere API query URL (see
#' \code{\link{get_url}}), keeps the timestamp and the requested variables
#' (rounded to two decimals), and writes the result to
#' \code{<path><id>.csv}. Called for its side effect; prints a message if
#' no data is available for the given station.
#' @param url query URL as returned by \code{\link{get_url}}
#' @param id station ID, used to construct the output file name
#' @param path output directory (character, with trailing slash)
#' @param var_names variable columns to keep (character vector; default \code{"rr"})
get_csv <- function(url, id, path, var_names = "rr")
{
  x <- tryCatch(read_csv(url, show_col_types = FALSE), error = function(e) NULL)
  if(is.null(x))
  {
    print(paste0("No data available for id ", id))
  }
  else
  {
    x <- x %>% dplyr::select(all_of(c("time", var_names)))
    x <- x %>% mutate(across(all_of(var_names), \(x) round(x, digits = 2)))
    path_out <- paste0(path, id, ".csv")
    print(path_out)
    write_csv(x, path_out)
  }
}

################################################################################
################################################################################