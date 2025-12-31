# Source: N/nwnr/Projekte/EROSA-Stat/Tawes_stations/download_station_data.R, 01.10.25
# split into script and functions (--> moved to R/dataPrep_download_functions.R)
################################################################################
######################## Download TAWES Stations ###############################
################################################################################
# inputs:
#   - API URL of the selected Geosphere Dataset (string)
#   - desired time frame to be covered (numeric vector, years)
#   - station IDs of selected stations (numeric vector)
# outputs:
#   - metadata on Austrain weather stations from Geosphere (csv)
#   - precipitation time series as csv
#       - column `time`: timestamp in 10min steps
#       - column `rr`: precipitation measured in each 10min step [mm]

library(tidyverse)
library(jsonlite)
devtools::load_all() #used functions: get_url, get_csv


### Get Metadata of all possible stations in Austria

url <- "https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min/metadata"
parsed_obj <- httr::content(httr::GET(url), as="text") %>% fromJSON()
stations <- parsed_obj$stations %>% as_tibble()
stations <- stations %>% mutate(valid_from = gsub("\\+00:00", "", valid_from), 
                                valid_to = gsub("\\+00:00", "", valid_to), valid_to = gsub("2100", "2022", valid_to))
stations <- stations %>% mutate(group_id = as.integer(group_id), id = as.integer(id))
climatic_variables <- parsed_obj$parameters[,1:4] # Precipitation variable name is rr
# We have to limit the stations dataset from 2009 to 2022 - this is the period where we have flash data
stations <- stations %>% filter(year(valid_from) <= 2009 & year(valid_to) == 2022) 

### Download data for the selected stations (Dornbirn=171; Bad Mitterndorf=10; Graz=30)
ids <- c(10,30,171)
stations <- stations %>% filter(id %in% ids)

write_csv(stations[c("id", "name", "state", "lon", "lat", "altitude")], "data/selected_stations_metadata.csv")

# Start of each time series can be set to 2009
stations <- stations %>% mutate(start = "2009-01-01T00:00:00")
# The url is added to the data frame so the stations can be downloaded
stations <- stations %>% mutate(url = pmap_chr(.l = list(id, start, valid_to), .f = function(id, start, valid_to)
  get_url(id = id, valid_from = start, valid_to = valid_to, parameter = "rr")))

# Set up directory where to save the data
path <- "data/raw_ts/"

if(!dir.exists(path))
{
  dir.create(path)
}

lapply(seq(length(ids)), function(ind) get_csv(url = stations$url[ind], id = stations$id[ind], path = path))

