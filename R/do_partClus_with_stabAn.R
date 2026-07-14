# Work horse function for partitioning clustering including global and segment
# level stability analysis

#' Function to run all partitioning clustering steps (repeated clustering over a
#' vector of k's; global stability analysis; and segment level stability analysis)
#' @return a list of:
#'          - global stability analysis results: ARI values (tibble with k and ARI)
#'          - segment level stabiity analysis results: ARI values (tibble with k, segment, ARI, and error_msg)
#'          - stepFlexclust object
#' @param df dataframe to be clustered
#' @param clusvars desired segmentation variables (character vector)
#' @param k numeric vector of k's to be tested
#' @param verbose Shall the functions give console status updates? Default=FALSE
#' @param b number of bootstrap runs
#' @examples
#' \dontrun{
#' # `scaled` as prepared in inst/scripts/004_partitioning_clustering.R
#' clusvars <- c('severity', 'magnitude', 'duration',
#'               'intensity', 'time_to_peak', 'flash')
#' part <- do_stabAn(scaled[["30"]], clusvars, k = 2:8, verbose = TRUE)
#' }
do_stabAn <- function(df, clusvars, k, verbose=FALSE, b=100) {

  res <- list(globalStability=NA,
              segLevelStability=NA,
              stepClustResults=NA)
  
  if(verbose) cat('Starting bootstraps!\n')
  
  boots <- flexclust::bootFlexclust(df[clusvars], k=k,
                                    nboot=b, verbose=verbose)
  res$globalStability <- boots@rand |> 
    as.data.frame() |> 
    pivot_longer(everything(),
                 names_to = 'k', values_to = 'ARI')
  
  if(verbose) cat('Starting stepclust!\n')
  
  res$stepClustResults <- flexclust::stepFlexclust(df[clusvars], k=k,
                                                   nrep=10, verbose=verbose) |> 
    flexclust::relabel()
  
  if(verbose) cat('Starting resamples!\n')
  
  resamples <- sapply(as.character(k), \(i) {
    
    if(verbose) cat(paste('k =', i, '\n'))
    
    tryCatch({
      flexclust::slswFlexclust(df[clusvars], res$stepClustResults[[i]])
    }, error = \(e) {
      conditionMessage(e)
    })
  }, simplify=FALSE)
  
  res$segLevelStability <- lapply(resamples, \(y) {
    if(data.class(y)=='resampleFlexclust') {
      t(y@validation[,1,]) |> 
        as.data.frame() |> 
        pivot_longer(everything(),
                     names_to = 'segment',
                     values_to = 'ARI') |> 
        mutate(error_msg=NA)
    } else {
      data.frame(
        segment='ERROR', ARI=NA,
        error_msg=as.character(y)
      )
    }
  }) |> 
    bind_rows(.id='k')
  
  res
  
}

if(FALSE) {
  do_stabAn(scaled$`30`, unname(clusvars), k=K, verbose=TRUE)
}