################################################################################
#### Aggregate theoretical distribution calculation & Cluster visualization ####
################################################################################
# inputs:
#  - filepath to precipitation event data objects for the desired stations in .RDS format created in step 3
#  - selected stations (numeric vector, length>=1)
#  - chosen Minimum Interevent Time (numeric vector, length=1)
#  - chosen filtering level (numeric vector, length=1)
#  - file path to clustering results in .RDA format created in step 5
#  - the chosen *k*s for these clustering results
#  - the (clustering) variables to be used in FAMD analysis (character vector)
#  - the desired dimensions for visualization in FAMD plot (list of 2 numeric vectors x and y)
#  - the covariates considered in PDP plots
#  - desired granularity of fitted theoretical distributions (numeric, length=1, of
#     desired output vector length)
# outputs:
#  - visualizations of the results:
#    - FAMD plots with cluster membership as color
#    - PDP plots for modelbased clustering
#    - ECDF plots (clusterwise and overall)
#    - plots of theoretical distribution fit (overall, and aggregate of clusterwise)
#  - aggregate theoretical distribution of precipitation, aggregated from clusterwise distributions
#    - comparison of distribution fits via Kullback-Leibler divergence

library(tidyverse)
library(flexmix)
library(lmom)
library(philentropy)
library(ggpubr)
library(FactoMineR)
devtools::load_all() #used functions: .get_clusters, aggregate_dist_calc, .get_KL, plt_famd, plt_dpd_flexmix, plt_ecdf, plt_cdf_comp

pth_prec <- .fix_path('') #TODO: replace this from script with hardcoded filepath in the next repo, data I guess
pth_res <- 'data/clusres'
iet <- 4
stations <- c(30, 171)
filt <- 0
k <- c('part30'=3,'part171'=3,
       'mod30'=2, 'mod171'=4)
clusvars <- c('severity', 'magnitude', 'duration',
              'intensity', 'time_to_peak', 'flash')
dims <- list(x=1, y=2:4)
covars <- c('magnitude', 'duration', 'time_to_peak')
length_out <- 1000

events <- paste0(pth_prec, 'model_input_events_iet_%dh_%d.RDS') |> 
  sprintf(iet, stations) |> 
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

clusses <- c(sprintf('part_iet_%d_station_%d.RDA',
                     iet, stations),
             sprintf('mod_iet_%d_filt_%g_station_%d.RDA',
                     iet, filt, stations))
for(i in clusses) paste0(pth_res, i) |> load()

clusses <- list(
  part30=list(dat=events$`30`,
              clus_short=.get_clusters(part_iet_4_station_30$full$stepClustResults,
                                       k['part30'])),
  part171=list(dat=events$`171`,
               clus_short=.get_clusters(part_iet_4_station_171$full$stepClustResults,
                                        k['part171'])),
  mod30=list(dat=events$`30`,
             mod=getModel(mod_iet_4_filt_0_station_30,
                          as.character(k['mod30'])),
             clus_short=.get_clusters(mod_iet_4_filt_0_station_30,
                                      k['mod30']),
             clus_long=.get_clusters(mod_iet_4_filt_0_station_30,
                                     k['mod30'],
                                     per_group=FALSE)),
  mod171=list(dat=events$`171`,
              mod=getModel(mod_iet_4_filt_0_station_171,
                           as.character(k['mod171'])),
              clus_short=.get_clusters(mod_iet_4_filt_0_station_171,
                                       k['mod171']),
              clus_long=.get_clusters(mod_iet_4_filt_0_station_171,
                                      k['mod171'],
                                      per_group=FALSE))
    
)

################################################################################
#---------------------Numeric comparison of distribution fit--------------------
################################################################################

KL <- sapply(clusses, \(obj) {
  rr <- obj$dat$rr
  if(is.null(obj$clus_long)) {
    obj$clus_long <- rep(obj$clus_short,
                         times=sapply(rr, \(x) sum(x>filt)))
  }
  pdf_theoretical <- 
    aggregate_dist_calc(rr, obj$clus_long,
                        filt=filt,
                        length_out = length_out,
                        type='PDF')
  rr <- unlist(rr)
  rr <- rr[rr>filt]
  pdf_empirical <- density(rr,
                           from = max(filt, 1e-3),
                           to = max(rr),
                           n = length_out)$y
  
  c(
    simple_Gamma=.get_KL(pdf_empirical, pdf_theoretical$simple_Gamma),
    aggregated_Gamma=.get_KL(pdf_empirical, pdf_theoretical$aggregated_Gamma)
  )
})

KL

################################################################################
#-----------------------Visualization of clustering results & ----------------------------
################################################################################


clusses$part30$lab <-'Graz Universität, Partitioning clustering, k=3'
clusses$part171$lab <-'Dornbirn, Partitioning clustering, k=3'
clusses$mod30$lab <-'Graz Universität, Clusterwise regression, k=2'
clusses$mod171$lab <- 'Dornbirn, Clusterwise regression, k=4'

## FAMD plots

lapply(clusses, \(obj) {
  plt_famd(obj$dat[clusvars], obj$clus_short,
           dims=dims, title=obj$lab)
}) %>% ggarrange(plotlist=.,
                 legend.grob=get_legend(.[[which.max(k)]] +
                                          theme(legend.position = 'bottom')),
                 legend='bottom')


## PDP plots
mod <- grep('mod', names(clusses))
lapply(clusses[mod], \(obj) {
  plt_dpd_flexmix(obj$mod, covars, fac_var='flash',
                  title=obj$lab)
}) %>% ggarrange(plotlist=., ncol=1,
                 legend.grob=get_legend(.[[which.max(k[mod])]] +
                                          theme(legend.position = 'bottom')),
                 legend='bottom')


## groupwise ECDFs
graz <- grepl('30', names(clusses))
foo <- function(x) {
  lbls <- clusses[x] |> lapply(`[[`, 'lab')
  station <- str_split_1(lbls[[1]], ',')[1]
  lbls <- setNames(str_split_fixed(lbls, ', ', 2)[,2],
                   names(lbls))
  list(station=station, lbls=lbls)
}
lbls_G <- foo(graz)
lbls_D <- foo(!graz)

list(
  plt_ecdf(events$`30`,
         clus_obj = lapply(clusses[graz], `[[`, 'clus_short'),
         title=lbls_G$station,
         subtitle=lbls_G$lbls),
  plt_ecdf(events$`171`,
         clus_obj = lapply(clusses[!graz], `[[`, 'clus_short'),
         title=lbls_D$station,
         subtitle=lbls_D$lbls)
) |> ggarrange(plotlist=_, ncol=1,
               common.legend=TRUE,
               legend='bottom')


## Visual comparison of distribution fit
list(
  plt_cdf_comp(events$`30`,
               clus_vec = lapply(clusses[graz], `[[`, 'clus_short'),
               filt=filt,
               title=lbls_G$station,
               subtitle=lbls_G$lbls),
  plt_cdf_comp(events$`171`,
             clus_vec = lapply(clusses[!graz], `[[`, 'clus_short'),
             filt=filt,
             title=lbls_D$station,
             subtitle=lbls_D$lbls)
) |> ggarrange(plotlist=_, ncol=1,
               common.legend=TRUE,
               legend='bottom')
