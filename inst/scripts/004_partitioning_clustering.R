################################################################################
##################### Partitioning clustering ##################################
################################################################################
# inputs:
#  - filepath to precipitation event data objects for the desired stations in .RDS format created in step 3
#  - selected stations (numeric vector, length>=1)
#  - chosen Minimum Interevent Time (numeric vector, length=1)
#  - starting seed
#  - segmentation variables to be used
#  - desired values of k to be tested
# outputs:
#  - RDA object containing clustering object, and ARIs from global and segment level stability analysis

library(tidyverse)
library(flexclust)
devtools::load_all() #used functions: do_stabAn

pth_prec <- 'data/events_combined/'
iet <- 4
stations <- c(171, 30)
start_seed <- as.numeric(as.Date('2025-09-01'))
K <- 2:8
clusvars <- c('severity', 'magnitude', 'duration',
              'intensity', 'time_to_peak', 'flash') |> 
  setNames(c('full', 'selected')[c(1,2,2,1,2,2)]) #TODO: this contains code for full and for selected. Show only the selected one in actually shared code


set.seed(start_seed)

events <- paste0(pth_prec,
                 'model_input_events_iet_', iet,
                 'h_*.RDS') |> 
  Sys.glob() %>% 
  grep(paste(stations, collapse='|'), x=_, value=TRUE) |> 
  sapply(readRDS, simplify=FALSE) %>%
  setNames(gsub('[^0-9]', '', names(.))) %>%
  setNames(gsub('4', '', names(.)))

scaled <- lapply(events, \(df) {
  mutate(df,
         flash = as.logical(flash),
         across(all_of(unname(clusvars[1:4])), log),
         across(all_of(unname(clusvars[1:5])), ~(. - min(.))/diff(range(.)))) |> 
    rename(precipitation=rr)
})

varnames <- sprintf('part_iet_%d_station_%d', iet, stations) |> 
  setNames(stations)

pth_out <- "data/clusres/"

if(!dir.exists(pth_out))
{
  dir.create(pth_out)
}

for(station in as.character(stations)) {
  part <- lapply(list(unname(clusvars),
                       unname(clusvars[names(clusvars)=='selected'])),
                  \(y) do_stabAn(scaled[[station]], y, k=K, verbose=TRUE)) |> 
    setNames(c('selected', 'full')) #TODO: This will simplify if only one set of variables is chosen
  assign(varnames[station], part)
  save(list=varnames[station],
       file=paste0(pth_out, varnames[station],
                   '.RDA'))
}

################################################################################
#-----------------Stability analysis: visualize results-------------------------
################################################################################
# all relevant outputs in this script have now been created, but let's look at the results of the stability analyses:

part_iet_4_station_30$full$stepClustResults
part_iet_4_station_30$full$globalStability |> 
  ggplot(aes(x=k, y=ARI)) + geom_boxplot()
part_iet_4_station_30$full$segLevelStability$error_msg |> 
  is.na() |> mean() #no errors, good
part_iet_4_station_30$full$segLevelStability |> 
  ggplot(aes(x=segment, y=ARI)) +
  geom_boxplot() + facet_wrap(~k, scales='free_x')
# and the same for the others...
