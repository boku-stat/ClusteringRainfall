# Source: N/nwnr/Projekte/EROSA-Stat/add_flash_and_precip.R. 01.10.25
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
devtools::load_all() #used functions: get_precip, get_flash_data

stations <- read_csv("data/selected_stations_metadata.csv")
ssf <- st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
# A buffer around the station of 7 km is used, in which each flash is counted
ssf_buffer <- st_transform(ssf, crs = 31287) %>% st_buffer(dist = units::set_units(7, "km"))

# Add flash data
flash <- read_csv("../ALDIS/ALDIS_csv/lightning.csv")
flash_sf <- st_as_sf(flash, coords = c("longitude", "latitude"), crs = 4326)
flash_sf <- st_transform(flash_sf, crs = 31287)
time <- tibble(date = seq(ymd_hms("2010-01-01 00:00:00"), ymd_hms("2022-12-31 23:00:00"), by = "hour"))
time <- time %>% mutate(month = month(date)) %>% filter(between(month, 4, 10)) %>% dplyr::select(all_of(c("date")))

joined <- st_join(ssf_buffer, flash_sf)
joined_sel <- joined %>% as_tibble() %>% dplyr::select(all_of(c("id","flash", "dt", "year", "month"))) %>% 
  filter(between(year, 2010, 2022) & between(month, 4, 10))

save_flash <- function(flash_data, id_sel, path_out)
{
  fname_out <- paste0(path_out, "flash_", id_sel, ".csv")
  xsel <- flash_data %>% filter(id == id_sel)
  xsel <- xsel %>% group_by(dt) %>% summarize(flash = ifelse(sum(flash) > 0, 1, 0))
  out <- xsel[c("flash", "dt")] %>% rename(date = dt)
  out <- left_join(time, out)
  out <- out %>% mutate(flash = ifelse(is.na(flash), 0, flash))
  write_csv(out, fname_out)
}

path_out <- "data/flash_events/"

if(!dir.exists(path_out))
{
  dir.create(path_out)
}

lapply(unique(joined_sel$id), save_flash, path_out = path_out, flash_data = joined_sel)

################################################################################
#--------- Joining Flash data, Event data and raw precipitation records --------
################################################################################

# List all event data
fl_events <- list.files("data/events/", pattern = "*.csv", full.names = TRUE)
fl_flash <- list.files("data/flash_events/", pattern = "*.csv", full.names = TRUE)
fl_raw <- list.files("data/raw_ts/", pattern = "*.csv", full.names = TRUE)

combine_data <- function(fname_event, path_out)
{
  bn <- basename(fname_event)
  id <- gsub(".csv", "", last(strsplit(bn, split = "_")[[1]]))
  ev <- read_csv(fname_event, show_col_types = FALSE)
  fname_raw <- grep(id, fl_raw, value = TRUE)
  x <- read_csv(fname_raw, show_col_types = FALSE)
  fname_flash <- grep(id, fl_flash, value = TRUE)
  flash <- read_csv(fname_flash, show_col_types = FALSE)
  # Add raw precipitation data to each event
  ev <- ev %>% mutate(rr = pmap(.l = list(start, end), .f = function(start, end)
    get_precip(start = start, end = end, df = x)))
  ev <- ev %>% mutate(flash = pmap_int(.l = list(start, end), .f = function(start, end) 
    get_flash_data(start = start, end = end, flash = flash)))
  fname_out <- paste0(path_out, "model_input_", bn)
  fname_out <- gsub(".csv", ".RDS", fname_out)
  saveRDS(ev, fname_out)
}

path_out <- "data/events_combined/"
if(!dir.exists(path_out))
{
  dir.create(path_out)
}

lapply(fl_events, combine_data, path_out = path_out)


