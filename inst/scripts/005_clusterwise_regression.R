################################################################################
########################## Clusterwise regression ##############################
################################################################################
# inputs:
#  - filepath to precipitation event data objects for the desired stations in .RDS format created in step 3
#  - selected stations (numeric vector, length>=1)
#  - chosen Minimum Interevent Time (numeric vector, length=1)
#  - starting seed
#  - model formula to be used
#  - desired values of k to be tested
#  - desired level of filtering of low precipitation values
#  - number of cores to run the bootstrapping analysis on
# outputs:
#  - RDA object containing clustering object, and ARIs from global and segment level stability analysis

library(tidyverse)
library(flexclust)
library(flexmix)
library(parallel)
devtools::load_all() #used functions: bootEventmix, modified MClapply (in script bootflexmix.R)

pth_prec <- .fix_path('') #TODO: replace this from script with hardcoded filepath in the next repo, data I guess
iet <- 4
stations <- c(30, 171)
start_seed <- as.numeric(as.Date('2025-08-26'))
K <- 2:8
clusvars <- c('severity', 'magnitude', 'duration',
              'intensity', 'time_to_peak', 'flash')
frm <- as.formula(precipitation + 10 ~ -1 + (log(magnitude+0.01) + log(duration+0.01) +
                                               time_to_peak)*flash)
filt <- c(0, 0.21)
cores <- max(1, floor(parallel::detectCores()*0.8))

set.seed(start_seed)

events <- paste0(pth_prec,
                 'model_input_events_iet_', iet,
                 'h_*.RDS') |> 
  Sys.glob() |> 
  grep(paste(stations, collapse='|'), x=_, value=TRUE) |> 
  sapply(readRDS, simplify=FALSE) %>%
  setNames(gsub('[^0-9]', '', names(.))) %>%
  setNames(gsub('4', '', names(.)))

scaled <- lapply(events, \(df) {
  mutate(df,
         flash = as.logical(flash),
         across(all_of(clusvars[1:4]), log),
         across(all_of(clusvars[1:5]), ~(. - min(.))/diff(range(.)))) |> 
    rename(precipitation=rr)
})

model_options <- expand.grid(list(iet=iet,
                             stationID=stations,
                             filt=filt))
varnames <- apply(model_options, 1,
                  \(x) sprintf('mod_iet_%d_filt_%g_station_%d',
                               x[1], x[3], x[2])) |> 
  gsub('.', '', x=_, fixed = TRUE)

pth_out <- "data/clusres/"

if(!dir.exists(pth_out))
{
  dir.create(pth_out)
}

for(i in 1:nrow(model_options)) {
  mod <- scaled[[as.character(model_options$stationID[i])]] |> 
    unnest(precipitation) |> 
    filter(precipitation > model_options$filt[i]) |> 
    stepFlexmix(y ~ x|id, k=K,
                data=_,
                model=FLXMRglm(frm, family='Gamma'),
                control=list(iter.max=100))
  assign(varnames[i], mod)
  save(list=varnames[i],
       file=paste0(pth_out, varnames[i],
                   '.RDA'))
}


################################################################################
#----------------------------Stability analysis---------------------------------
################################################################################

varnames <- paste(varnames, 'stabAn', sep='_')

#--DON'T RUN----------------
# this section was run on a HPC, using roughly 30 cores, taking about a week.
# Run this full setup only on an appropriate computing setup.
for(i in 1:nrow(model_options)) {
  mod <- scaled[[as.character(model_options$stationID[i])]] |> 
    bootEventmix(k=K, formula=y~x|id,
                 filter=model_options$filt[i],
                 model=FLXMRglm(frm, family='Gamma'),
                 nboot=100,
                 control=list(iter.max=100),
                 multicore=TRUE,
                 mc.cores=cores)
  assign(varnames[i], mod)
  save(list=varnames[i],
       file=paste0(pth_out, varnames[i],
                   '.RDA'))
}

################################################################################
#-----------------Stability analysis: visualize results-------------------------
################################################################################
# all relevant outputs in this script have now been created, but let's look at the results of the global stability analysis
mod <- list(db_0=mod_iet_4_filt_0_station_171_stabAn$rand,
            db_021=mod_iet_4_filt_021_station_171_stabAn$rand,
            g_0=mod_iet_4_filt_0_station_30_stabAn$rand,
            g_021=mod_iet_4_filt_021_station_30_stabAn$rand) |> 
  lapply(as.data.frame) |> 
  lapply(pivot_longer, everything(),
         names_to='k', values_to='ARI') |> 
  bind_rows(.id='stat_filt') |> 
  mutate(station=str_split_fixed(stat_filt, '_', 2)[,1],
         filtering=ifelse(grepl("21", stat_filt), 0.21, 0)) |> 
  select(-stat_filt)

ct <- mutate(mod, ran=!is.na(ARI)) |> group_by(k, station, filtering) |>
  reframe(successes=sum(ran), runs=n()) |> 
  mutate(successes_text=paste0(successes, '%'),
         ARI=0.22)

left_join(mod, select(ct, -c(ARI, runs, successes_text))) |> 
  ggplot(aes(x=k, y=ARI)) +
  geom_boxplot(aes(fill=successes)) +
  geom_text(data=ct, aes(label=successes_text),
            col=formals(scale_fill_gradient)$low) +
  facet_grid(station~filtering)
