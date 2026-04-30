################################################################################
######################## MASTER ANALYSIS SCRIPT ################################
################################################################################
# This script orchestrates the complete analysis pipeline for multiple minimum
# interevent time (MIT) values.Specifically, it documents our robustness analyses on
# interevent time values different to our baseline value of 4hrs.
# We run our robustness analysis on MIT (also called "IET" in some parts of the code)
# values of 3, 6, 12, 18 and 24hrs.
#
# The script runs, for each MIT value:
#   1. Precipitation data download (note that flash data is behind a paywall and must
#      thus be procured and downloaded in a separate application step.)
#   2. Data preparation (flash + precipitation)
#   3. Partitioning clustering with stability analysis
#   4. Clusterwise regression with stability analysis
#   5. Report generation (graphs and Kullback-Leibler divergence values)
#
################################################################################

library(tidyverse)
library(jsonlite)
library(zoo)
library(sf)
library(tidync)
library(flexclust)
library(flexmix)
library(parallel)
library(rmarkdown)
library(philentropy)
devtools::load_all() #alternatively, source(...) #TODO: add relevant function script paths here

################################################################################
# SETUP
################################################################################

# IET values for robustness analysis
IET_VALUES <- c(3, 6, 12, 18, 24)
# Station selection
STATIONS <- c(30, 171)
# Time series range
RANGE <- list()
RANGE$START <- ymd_hms("2010-01-01 00:00:00")
RANGE$END <- ymd_hms("2022-12-31 23:50:00")
RANGE$MONTHS <- c(4, 10)
# Event start thresholds
THRESHOLDS <- list()
THRESHOLDS$INIT <- 0.2
THRESHOLDS$SUM <- 1.27
# Lightning radius
FLASH_RADIUS <- units::set_units(7, "km")
# Seeds for reproducibility
SEED <- as.numeric(as.Date('2026-01-08'))
# Clustering parameters
K_VALUES <- 2:8
CLUS_VARS <- c('severity', 'magnitude', 'duration',
               'intensity', 'time_to_peak', 'flash')
# Clusterwise regression parameters
REGRESSION_FORMULA <- as.formula(
  precipitation + 10 ~ -1 + (log(magnitude+0.01) + log(duration+0.01) +
                              time_to_peak)*flash
)
FILTER_VALUES <- 0 #TODO (note to self): Only ran for filt0
# Computational settings
N_CORES <- max(1, floor(parallel::detectCores() * 0.8))
N_BOOTSTRAP <- 100  # reduce for faster testing, 100+ for production
VERBOSE <- TRUE
MClapply_wrapper <- function(X, FUN, ..., multicore=TRUE) {
  ClusteringRainfall:::MClapply(X, FUN, ...,
                                multicore=multicore,
                                mc.cores=N_CORES)
}
assignInNamespace("MClapply", MClapply_wrapper, ns="flexclust") ##overriding flexclust's default of 2 cores
# File paths
PATHS <- list()
PATHS$STATION_METADATA <- "data/selected_stations_metadata.csv"
PATHS$FLASH_RAW <- "data/ALDIS/lightning.csv" #note that this data is behind a paywall and must be added separately
PATHS$EVENTS <- "data/events/"
PATHS$FLASH_EVENTS <- "data/flash_events/"
PATHS$RAW_TS <- "data/raw_ts/"
PATHS$EVENTS_COMBINED <- "data/events_combined/"
PATHS$RESULTS <- "data/clusres/"
PATHS$REPORTS <- "inst/reports/"
#TODO: check whether I need all of those
# Create directories if they don't exist
for (path in PATHS) {
  if (grepl("/$", path) & !dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

#******
#hiiiiiii - for the reports, change afterwards
IET_VALUES <- 3
PATHS$EVENTS_COMBINED <- '/media/nwnr/Projekte/EROSA-Stat/Tawes_stations/data/events_combined/'
PATHS$RESULTS <- '/media/nwnr/Projekte/EROSA-Stat/Tawes_stations/data/clusres/'
CHOSEN_Ks <- c(3, 3, 2, 2) |> 
  setNames(c('part171', 'part30', 'mod171', 'mod30'))
#

################################################################################
# STEP 1: DOWNLOAD PRECIPITATION TIME SERIES DATA
################################################################################
# completes steps from script 001_download_station_data.R

step1_download_station_data <- function(ids) {
  if(VERBOSE) message(sprintf("\n=== STEP 1: Download station data for stations %s===\n",
                              paste(ids, collapse=", ")))
  
  ### Get Metadata of all possible stations in Austria
  
  url <- "https://dataset.api.hub.geosphere.at/v1/station/historical/klima-v2-10min/metadata"
  parsed_obj <- httr::content(httr::GET(url), as="text") %>% fromJSON()
  stations <- parsed_obj$stations %>% as_tibble()
  stations <- stations %>% mutate(valid_from = gsub("\\+00:00", "", valid_from), 
                                  valid_to = gsub("\\+00:00", "", valid_to), valid_to = gsub("2100", "2022", valid_to))
  stations <- stations %>% mutate(group_id = as.integer(group_id), id = as.integer(id))
  climatic_variables <- parsed_obj$parameters[,1:4] # Precipitation variable name is rr
  stations <- stations %>% filter(year(valid_from) <= year(RANGE$START) & 
                                    year(valid_to) == year(RANGE$END)) 
  
  ### Download data for the selected stations (Dornbirn=171; Bad Mitterndorf=10; Graz=30)
  stations <- stations %>% filter(id %in% ids)
  
  write_csv(stations[c("id", "name", "state", "lon", "lat", "altitude")], "data/selected_stations_metadata.csv")
  
  # Start of each time series can be set to 2009
  stations <- stations %>% mutate(start = "2009-01-01T00:00:00")
  # The url is added to the data frame so the stations can be downloaded
  stations <- stations %>% mutate(url = pmap_chr(.l = list(id, start, valid_to), .f = function(id, start, valid_to)
    get_url(id = id, valid_from = start, valid_to = valid_to, parameter = "rr")))

  # Download precipitation data for the selected stations
  lapply(seq(length(ids)), function(ind) {
    get_csv(url = stations$url[ind],
            id = stations$id[ind],
            path = PATHS$RAW_TS)
  })

}

# Usage: step1_download_station_data(ids=STATIONS)

################################################################################
# STEP 2: SEPARATE PRECIPITATION TIME SERIES INTO EVENTS
################################################################################
# completes steps from script 002_create_event_data.R

step2_create_event_data <- function(iets) {
  if(VERBOSE) message(sprintf("\n=== STEP 2: Separate Precipitation Time Series Into Events (IETs=%s) ===\n",
                              paste(iets, collapse=", ")))
  fl <- list.files(path = PATHS$RAW_TS,
                   pattern = "*.csv", full.names = TRUE)
  df <- tibble(date = seq(RANGE$START, RANGE$END, by = "10 mins"))
  
  for(iet in iets) {
    lapply(fl, get_events, 
           df=df,nh=iet,
           path_out=PATHS$EVENTS,
           th_init=THRESHOLDS$INIT,
           th_sum=THRESHOLDS$SUM)    
  }

}

# Usage: step2_create_event_data(iets=IET_VALUES)

################################################################################
# STEP 3: ADD LIGHTNING OCCURRENCE TO EVENTS
################################################################################
# completes steps from script 003_add_flash_and_precip.R

step3_add_flash_and_precip <- function() {
  if(VERBOSE) message("\n=== STEP 3: Add Lightning Occurrence To Events ===\n")
  
  stations <- read_csv(PATHS$STATION_METADATA)
  ssf <- st_as_sf(stations, coords = c("lon", "lat"), crs = 4326)
  ssf_buffer <- st_transform(ssf, crs = 31287)  |> 
    st_buffer(dist = FLASH_RADIUS)
  
  # Add flash data
  flash <- read_csv(PATHS$FLASH_RAW)
  flash_sf <- st_as_sf(flash, coords = c("longitude", "latitude"), crs = 4326)
  flash_sf <- st_transform(flash_sf, crs = 31287)
  time <- tibble(date = seq(RANGE$START, RANGE$END, by = "hour"))
  time <- time |> mutate(month = month(date)) |>
    filter(between(month, RANGE$MONTHS[1],
                   RANGE$MONTHS[2])) |> dplyr::select(all_of(c("date")))
  
  joined <- st_join(ssf_buffer, flash_sf)
  joined_sel <- joined |> as_tibble() |> dplyr::select(all_of(c("id","flash", "dt", "year", "month"))) |> 
    filter(between(year, year(RANGE$START), year(RANGE$END)) 
           & between(month, RANGE$MONTHS[1],
                     RANGE$MONTHS[2]))
  
  lapply(unique(joined_sel$id), save_flash,
         path_out = PATHS$FLASH_EVENTS,
         flash_data = joined_sel,
         time_data = time)

  # combine data sources
  fl_events <- list.files(PATHS$EVENTS, pattern = "*.csv", full.names = TRUE)
  fl_flash <- list.files(PATHS$FLASH_EVENTS, pattern = "*.csv", full.names = TRUE)
  fl_raw <- list.files(PATHS$RAW_TS, pattern = "*.csv", full.names = TRUE)
  
  lapply(fl_events, combine_data,
         path_out = PATHS$EVENTS_COMBINED,
         fname_flash = fl_flash,
         fname_precip = fl_raw)
}

# Usage: step3_add_flash_and_precip()


################################################################################
# STEP 4: PARTITIONING CLUSTERING
################################################################################
# completes steps from script 004_partitioning_clustering.R

step4_partitioning_clustering <- function(iet) {
  if(VERBOSE) message(sprintf("\n=== STEP 4: Partitioning Clustering (IET=%d) ===\n", iet))
  
  set.seed(SEED)

  # Load event data
  events <- sprintf("%smodel_input_events_iet_%dh_*.RDS",
                    PATHS$EVENTS_COMBINED, iet) |> 
    Sys.glob() %>% 
    grep(paste(STATIONS, collapse='|'), x = ., value = TRUE) |> 
    sapply(readRDS, simplify = FALSE) %>%
    setNames(gsub('.*_(\\d+)\\.RDS$', '\\1', names(.)))
  
  if (length(events) == 0) {
    warning(sprintf("No event data found for IET=%d", iet))
    return(NULL)
  }
  
  # Scale variables
  scaled <- lapply(events, function(df) {
    mutate(df,
           flash = as.logical(flash),
           across(all_of(CLUS_VARS[1:4]), log),
           across(all_of(CLUS_VARS[1:5]), ~(. - min(.)) / diff(range(.)))) %>%
      rename(precipitation = rr)
  })
  
  # Run clustering for each station
  for (stn in as.character(STATIONS)) {
    
    varname <- sprintf('part_iet_%d_station_%s', iet, stn)
    
    part <- do_stabAn(scaled[[stn]], CLUS_VARS, k = K_VALUES,
                      verbose = VERBOSE,
                      b = N_BOOTSTRAP)
    
    fname_out <- sprintf("%s%s.RDA", PATHS$RESULTS, varname)
    save(list = "part", file = fname_out)
  }
  
}

# Usage:  step4_partitioning_clustering(iet=max(IET_VALUES))


################################################################################
# STEP 5: CLUSTERWISE REGRESSION
################################################################################
# completes steps from script 005_clusterwise_regression.R

step5_clusterwise_regression <- function(iet) {

  if(VERBOSE) message(sprintf("\n=== STEP 5: Clusterwise Regression (IET=%d) ===\n", iet))
  
  set.seed(SEED)
  
  # Load event data
  events <- sprintf("%smodel_input_events_iet_%dh_*.RDS",
                    PATHS$EVENTS_COMBINED, iet) |> 
    Sys.glob() %>% 
    grep(paste(STATIONS, collapse='|'), x = ., value = TRUE) |> 
    sapply(readRDS, simplify = FALSE) %>%
    setNames(gsub('.*_(\\d+)\\.RDS$', '\\1', names(.)))
  
  if (length(events) == 0) {
    warning(sprintf("No event data found for IET=%d", iet))
    return(NULL)
  }
  
  # Scale variables
  scaled <- lapply(events, function(df) {
    mutate(df,
           flash = as.logical(flash),
           across(all_of(CLUS_VARS[1:4]), log),
           across(all_of(CLUS_VARS[1:5]), ~(. - min(.)) / diff(range(.)))) %>%
      rename(precipitation = rr)
  })
  
  # Create model options grid
  model_options <- expand.grid(list(iet = iet,
                                   stationID = STATIONS,
                                   filt = FILTER_VALUES))
  
  varnames <- apply(model_options, 1,
                   function(x) sprintf('mod_iet_%d_filt_%g_station_%d',
                                     x[1], x[3], x[2])) %>%
    gsub('.', '', x = ., fixed = TRUE)
  
  # Run regression models
  for (i in 1:nrow(model_options)) {
    station_id <- model_options$stationID[i] |> 
      as.character()
    filt_val <- model_options$filt[i]
    varname <- varnames[i]
    
    mod <- scaled[[station_id]] |> 
      unnest(precipitation) |> 
      filter(precipitation > filt_val) |> 
      stepFlexmix(y ~ x|id, k = K_VALUES,
                 data = _,
                 model = FLXMRglm(REGRESSION_FORMULA, family = 'Gamma'),
                 control = list(iter.max = 100))
    
    fname_out <- sprintf("%s%s.RDA", PATHS$RESULTS, varname)
    save(list = "mod", file = fname_out)
    
    # Bootstrap stability analysis
    stabAn <- scaled[[station_id]] |> 
      bootEventmix(k = K_VALUES, 
                    formula = y ~ x|id,
                    filter = filt_val,
                    model = FLXMRglm(REGRESSION_FORMULA, family = 'Gamma'),
                    nboot = N_BOOTSTRAP,
                    control = list(iter.max = 100),
                    multicore = TRUE,
                    mc.cores = N_CORES)
      
    stabname <- paste0(varnames[i], '_stabAn')
    results[[stabname]] <- stabAn
      
    fname_out <- sprintf("%s%s.RDA", PATHS$RESULTS, stabname)
    save(list = "stabAn", file = fname_out)
    
  }
}

# Usage: step5_clusterwise_regression(iet=max(IET_VALUES))

################################################################################
# STEP 6: AGGREGATE DISTRIBUTION FIT, AND REPORT GENERATION
################################################################################
# completes steps from script 006_aggregate_distribution_calc.R, and collects
# the visual and numeric comparisons in a report

step6_generate_report <- function(iet) {
  if(VERBOSE) message(sprintf("\n=== STEP 6: Generating Report (IET=%d) ===\n", iet))
  
  report_file <- "inst/scripts/101_analysis_report_template.Rmd"
  output_file <- sprintf("report_iet_%d.pdf", iet)
  
  regvars <- REGRESSION_FORMULA |> 
    as.character() |> 
    sapply(CLUS_VARS, grepl, x=_) |> 
    colSums() |> as.logical()
  
  tryCatch({
    rmarkdown::render(
      input = report_file,
      output_file = output_file,
      output_dir = PATHS$REPORTS,
      params = list(iet = iet,
                   stations = STATIONS,
                   path_results = PATHS$RESULTS,
                   path_dat = PATHS$EVENTS_COMBINED,
                   filt=FILTER_VALUES,
                   ks=K_VALUES,
                   clusvars=CLUS_VARS,
                   regvars=CLUS_VARS[regvars],
                   chosen_ks=CHOSEN_Ks), #default setting: 'NULL'
      envir = new.env()
    )
   
  }, error = function(e) {
    warning(sprintf("Error generating report for IET=%d: %s", iet, e$message))
  })
}

#Usage: step6_generate_report(iet=4) #still not yet fully tested

################################################################################
# MAIN EXECUTION
################################################################################

main <- function(step1=TRUE,
                 step2=TRUE,
                 step3=TRUE,
                 step4=TRUE,
                 step5=TRUE,
                 step6=TRUE) {
  
  if (step1) {
    step1_download_station_data(ids=STATIONS)
  } else {
    message("\n=== STEP 1: Skipped ===\n")
  }
  
  if (step2) {
    step2_create_event_data(iets=IET_VALUES)
  } else {
    message("\n=== STEP 2: Skipped ===\n")
  }
  
  if (step3) {
    step3_add_flash_and_precip()
  } else {
    message("\n=== STEP 3: Skipped ===\n")
  }
  
  for (iet in IET_VALUES) {
    
    if (step4) {
      step4_partitioning_clustering(iet=iet)
    } else {
      message("\n=== STEP 4: Skipped ===\n")
    }
    
    if (step5) {
      step5_clusterwise_regression(iet=iet)
    } else {
      message("\n=== STEP 5: Skipped ===\n")
    }
    
    if (step6) {
      step6_generate_report(iet=iet)
    } else {
      message("\n=== STEP 6: Skipped ===\n")
    }
    
    
    
  }
  
}

################################################################################
# RUN THE PIPELINE
################################################################################
# Run the selected pipeline steps

main(
  step1 = FALSE, # download precipitation time series data
  step2 = FALSE, # derive events and event characteristics
  step3 = FALSE, # combine event and flash data
  step4 = FALSE, # partitioning clustering (+stability analysis)
  step5 = FALSE, # clusterwise regression (+stability analysis)
  step6 = TRUE # generate diagnostic reports
)

