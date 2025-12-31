# Source: N/nwnr/Projekte/EROSA-Stat/Tawes_stations/create_event_data.R, 01.10.25
################################################################################
##################### Get events for 10-min precipitation events ###############
################################################################################
# inputs:
#   - filepath to the precipitation time series .csvs downloaded in step 1
#   - desired Minimum Inter Event Time values (numeric vector)
# outputs: list of events for each evaluated station, plus event characteristics (yet without lightning) in .csv
#   - columns: event ID, start & end time [time, 10min timepoint], severity [mm], magnitude [mm], duration [min], intensity [mm/min], time_to_peak


library(tidyverse)
library(zoo)
devtools::load_all() #used functions: get_events (calls on get_event_characteristics, and get_max_hourly_precip)

################################################################################

fl <- list.files(path = "data/raw_ts", pattern = "*.csv", full.names = TRUE)
start <- ymd_hms("2010-01-01 00:00:00")
end <- ymd_hms("2022-12-31 23:50:00")
df <- tibble(date = seq(start, end, by = "10 mins"))

# Generate Event characteristics - Using different Interevent time criteria

path_out <- "data/events/"

if(!dir.exists(path_out))
{
  dir.create(path_out)
}

iet <- c(3,4,6,12,18,24)
for(i in iet)
{
  lapply(fl, get_events, path_out = path_out, th_init = 0.2, nh = i, 
         th_sum = 1.27)
}