
# Source: N/nwnr/Projekte/EROSA-Stat/Tawes_stations/download_station_data.R, 01.10.25
# split into functions and script (<-- in inst/001_download_station_data.R)


################################################################################
########## Functions for downloading TAWES stations data #######################
################################################################################

get_url <- function(id, valid_from, valid_to, parameter)
{
  url <- paste0("https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min?parameters=", parameter, "&start=", valid_from, 
                "&end=", valid_to, "&station_ids=", id, "&output_format=csv")
  url
}

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