#Functions for partitioning clustering
#Stand: 05.02.25

#' Function to fit&choose the pam model for whose `k` avg. silhouette width is maximal
#' @return the `pam` object for whose `k` the avg. silhouette width is maximal
#' @param distmat a distance matrix. (Would also work on x directy)
#' @param K a numeric vector declaring the number of clusters to test (minimum
#'          value of 2)
silchoice <- function(distmat, K) {
  
  res <- setNames(lapply(K, \(k) cluster::pam(x=distmat, k=k)), K)
  bestk <- sapply(res, \(y) y$silinfo$avg.width) |> 
    which.max() |> names()
  res[[bestk]]
  
}

#' Function to map each clustering to the desired descriptive variables from precip
#' @return a list of dataframes of `length=length(stations)` containing the desired
#'         descriptive variables, plus the clustering
#' @param cluslist a list of clustering models of length `unique(dat$station)`.
#'        Either of class `pam` or `ClusterMixedDataModel`
#' @param dat the data object with which to merge. Provide either a path to a csv
#'        file, or a data frame which is already read in. Default='data/precipitationEvents.csv'.
#'        Needs to have the column 'station'.
#' @param cols a character vector of the columns in dat with which to merge.
#'        Default=`dplyr::everything()`
#' @param stationnames the desired names of the output list items. Default=NULL
#'        (<-- defaults to an unnamed list)
#' @param relabelBy variable by which to relabel the clusters. After all, labelling
#'        is arbitrary, and ordering of clusters in plots will not be consistent.
#'        Default=NULL: no relabelling occurs.
extract_clusters <- function(cluslist, dat='data/precipitationEvents.csv',
                             cols=dplyr::everything(),
                             stationnames=NULL,
                             relabelBy=NULL) {
  cl <- sapply(cluslist, data.class)
  
  if(all(cl=='pam')) {
    extr <- function(mod) mod$clustering
  } else if(all(cl=='ClusterMixedDataModel')) {
    extr <- function(mod) slot(mod, 'zi')
  } else {stop('Oh oh!')}
  
  dat <- .split_precip(file=dat, cols=cols, stationnames=stationnames)
  clus <- sapply(cluslist, extr, simplify = F) |> 
    sapply(as.factor)
  
  if(!is.null(relabelBy)) {
    clus <- Map(\(x, y) {
      fct_reorder(y, x[[relabelBy]], .fun=mean) |> #, .desc=TRUE) |>
        fct_relabel(~as.character(seq_along(.)))
      }, x=dat, y=clus)
  }
  
  Map(\(x, y) cbind(x, cluster=y), x=dat, y=clus)
   
}

#' Get the scaled mean value for each station (clusterwise and total) for desired variables
#' @return a data.frame where each row shows the mean values by region and cluster,
#'         and total mean values for each region
#' @param mergedat a dataframe to which the 'cluster' column has already been merged.
#' @param cols the columns for which means shall be calculated. Numeric columns
#'        are scaled previously
#' @param group_vars the grouping variables, one of which should be 'cluster'.
#' @param longer shall the result be converted to a longer format? Default=TRUE
get_clusterMeans <- function(mergedat,
                             cols=c('severity', 'magnitude',
                                    'duration', 'intensity',
                                    'time2peak', 'flash'),
                             group_vars=c('mainregion', 'cluster'),
                             longer=T) {
  
  get_Means <- function(dat=mergedat, grpvrs=group_vars,
                        total=c(T,F)) {
    if(total) {
      grpvrs <- setdiff(group_vars, 'cluster')
    } else { grpvrs <- group_vars }
    
    select(dat, all_of(c(grpvrs, cols))) |> 
      mutate(across(where(is.numeric), scale)) |> 
      group_by(across(all_of(grpvrs))) |> 
      reframe(across(cols, mean))
  }
  
  res <- get_Means(total=F) |> 
    left_join(get_Means(total=T),
              by=setdiff(group_vars, 'cluster'),
              suffix=c('', '.total'))
  
  if(longer) res <- pivot_longer(res, where(is.numeric)) |> 
    mutate(type=ifelse(endsWith(name, 'total'),
                       'total', 'clusterwise'),
           name=gsub('.total', '', name)) |> 
    pivot_wider(id_cols=c(all_of(group_vars), name),
                names_from = type)
  
  res
}

#' Takes factor (or logical) variables from a grouped data set, returns proportional
#' counts for each factor level, and widens it to columns
#' @param data the data set. Needs to contain at least one grouping variable, and
#'        the factor to count
#' @param group_vars a character vector of variables by which to group (obviously
#'        need to be part of `data`)
#' @param count_var the factor/logical variable in `data` where the level proportions
#'        shall be counted
summarise_proportions <- function(data, group_vars, count_var) {
  res <- data |> 
    group_by(across(all_of(group_vars))) |> 
    count(.data[[count_var]]) |> 
    pivot_wider(names_from = all_of(count_var), values_from = n, values_fill = 0) |> 
    ungroup()
  cls <- setdiff(colnames(res), group_vars)
  mutate(res, tot=rowSums(res[,cls]),
         across(all_of(cls), ~.x/tot)) |> 
    select(-tot)
}

#' Get model names of MixAll::clusterMixedDataModel and paste by '...'
#' (same structure I used to set the names), plus ordered with gamma first
#' #basically summary() without the extra info
#' @param mixmod_object object of class MixAll::clusterMixedDataModel
get_modnames <- function(mixmod_object) {
  lapply(mixmod_object@lcomponent, slot, 'modelName') |> 
    rev() |> paste(collapse='...')
}

#' Get model specifications available in MixAll
#' @param fun function from MixAll::...Names
# Find all object assignments (lists of model names)
extract_MixAll_modspecs <- function(fun) {

  fun_nm <- deparse(substitute(fun))
  if(grepl('Gamma', fun_nm)) {
    modname <- 'gamma'
  }
  else if(grepl('Categorical', fun_nm)) {
    modname <- 'categorical'
  } 
  else stop('Error, only implemented for clusterGammaNames and clusterCategoricalNames')
  
  bod <- as.list(body(fun))
  assignments <- sapply(bod, \(y) any(grepl(modname, y)))
  bod <- bod[assignments]
  
  vals <- lapply(bod, eval)
  nams <- sapply(bod, \(y) as.character(y[[2]]))
  names(vals) <- nams
  
  return(vals)
}
deparse(substitute(MixAll::clusterGammaNames))
#fix bootstrap thing when MixAll is running

#' Extract clusters from flex model with desired k
#' This will work both on stepFlexclust and stepFlexmix objects;
#' and, for flexmix objects, 
#' @param clus_obj clustering object of class `stepFlexclust` or `stepFlexmix`
#' @param k desired number of clusters
#' @param per_group Only used for `stepFlexmix` objects: Shall the function return one
#'                  element of the clustering vector per group, or shall it return the whole
#'                  clustering vector? (This param is irrelevant in case of no grouping, as
#'                  groups will simply be a vector of unique levels.)
.get_clusters <- function(clus_obj, k, per_group=TRUE) {
  if(is.numeric(k)) k <- as.character(k)
  mod <- modeltools::getModel(clus_obj, k)
  clus <- modeltools::clusters(mod)
  if(data.class(clus_obj)=="stepFlexmix" & per_group) {
    grp <- flexmix::group(mod)
    clus <- clus[!duplicated(grp)] #irrelevant whether grouping was used or not, groups will just be unique in case of no grouping
  } #side note: this ofc assumes that the clusters for each group will be equal. This applies in our case, I'm unsure whether there are applications where that is not the case (would require to add error testing for that case)
  clus
}
 