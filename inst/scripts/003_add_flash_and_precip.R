################################################################################
############ Add flash and precipitation data to events#########################
################################################################################
# inputs:
#   - filepaths to: station metadata file (csv); raw lightning data (csv);
#                   raw precipitation data downloaded in step 1 (csv);
#                   event level data created in step 2 (csv).
#   - desired time frame to be covered (vector of length=2 with start and end date as datetime)
# outputs:  list of events for each evaluated station, plus event characteristics including lightning and also nested raw precipitation time series in .RDS
#   - columns: event ID, start & end time [time, 10min timepoint], severity [mm], magnitude [mm], duration [min], intensity [mm/min], time_to_peak,
#              flash [TRUE/FALSE]; rr [precipitation, nested list of mm of precipitation measured in each 10min timestep]


library(tidyverse)
library(sf)
library(tidync)
devtools::load_all() #used functions: get_precip, get_flash_data, save_flash, combine_data

stations <- read_csv("data/selected_stations_metadata.csv")
ssf <- st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
# A buffer around the station of 7 km is used, in which each flash is counted
ssf_buffer <- st_transform(ssf, crs = 31287)  |> st_buffer(dist = units::set_units(7, "km"))

# Add flash data
#--DON'T RUN--- caution: raw file needed for the next line is not provided as it is owned by ALDIS. Processed stationwise flash data (as produced by this script, see below) is provided in this repository
flash <- read_csv("data/ALDIS/lightning.csv")
#-------
flash_sf <- st_as_sf(flash, coords = c("longitude", "latitude"), crs = 4326)
flash_sf <- st_transform(flash_sf, crs = 31287)
time <- tibble(date = seq(ymd_hms("2010-01-01 00:00:00"), ymd_hms("2022-12-31 23:00:00"), by = "hour"))
time <- time |> mutate(month = month(date)) |> filter(between(month, 4, 10)) |> dplyr::select(all_of(c("date")))

joined <- st_join(ssf_buffer, flash_sf)
joined_sel <- joined |> as_tibble() |> dplyr::select(all_of(c("id","flash", "dt", "year", "month"))) |> 
  filter(between(year, 2010, 2022) & between(month, 4, 10))

path_out <- "data/flash_events/"

if(!dir.exists(path_out))
{
  dir.create(path_out)
}

lapply(unique(joined_sel$id), save_flash, 
       path_out = path_out,
       flash_data = joined_sel,
       time_data = time)

################################################################################
#--------- Joining Flash data, Event data and raw precipitation records --------
################################################################################

# List all event data
fl_events <- list.files("data/events/", pattern = "*.csv", full.names = TRUE)
fl_flash <- list.files("data/flash_events/", pattern = "*.csv", full.names = TRUE)
fl_raw <- list.files("data/raw_ts/", pattern = "*.csv", full.names = TRUE)

path_out <- "data/events_combined/"
if(!dir.exists(path_out))
{
  dir.create(path_out)
}

lapply(fl_events, combine_data,
       path_out = path_out,
       fname_flash=fl_flash,
       fname_precip=fl_raw)

#caution, moved the functions, check
