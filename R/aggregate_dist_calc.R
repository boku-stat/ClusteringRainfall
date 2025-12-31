#' Helper function for fitting a Gamma distribution to a vector of precipitation values
#' @param prec_actual precipitation vector from data (filtered to values bigger 0 only)
#' @param prec_predict the 'newdata' precipitation vector, on which the Gamma distribution
#'                     fitted onto `prec_actual` shall be predicted
#' @param type Character value of either `PDF` or `ECDF`. Shall the probability density function
#'             of `prec_predict` be predicted, or the empirical cummulative distribution function?
.calc_gamma <- function(prec_actual, prec_predict,
                        type=c('PDF', 'ECDF')) {
  
  params <- lmom::samlmu(prec_actual) |> 
    lmom::pelgam()
  if(type=='ECDF') lmom::cdfgam(prec_predict, para=params)
  else if(type=='PDF') dgamma(prec_predict,
                              shape=params['alpha'],
                              scale=params['beta']) #see: ?cdfgam
  else (stop('Wrong type parameter!'))
  
}

#' Function to fit Gamma distribution to precipitation data - both over all data points,
#' and by cluster (for `stepFlexclust` and `stepFlexmix` objects).
#' @param precip_list a list of precipitation vectors, nested by precipitation event (Caution:
#'                    in case of clusterwise regression, this shall be the same input as entered
#'                    into the clusterwise regression.)
#' @param filt the filtering value for the precipitation values previous to Gamma fitting (Caution:
#'                    in case of clusterwise regression, this must be the same input as entered
#'                    into the clusterwise regression.)
#' @param clus_obj If the clustering vector hasn't been extracted yet:
#'                 a clustering object of class `stepFlexclust` or `stepFlexmix`.
#'                 If it has: vector of cluster assignments.
#' @param k the desired number of clusters out of the `k` values tested in the clustering object
#'          (not used if `clus_obj` is the already extracted clustering vector)
#' @param length_out length of prediction vector. The limits of the prediction vector are determined by
#'                   the values in the `precip_list` and/or `filt`, this will only control the "smoothness of the predicted curve"
#' @param type Character value of either `PDF` or `ECDF`. Shall the probability density function
#'             of `prec_predict` be predicted, or the empirical cummulative distribution function?
aggregate_dist_calc <- function(precip_list, filt,
                                clus_obj, k=NULL,
                                length_out=1000,
                                type=c('PDF', 'ECDF')) {

  if(data.class(clus_obj) %in% c('stepFlexclust', 'stepFlexmix')) {
    if(is.null(k)) stop('Please specify desired number of clusters k!')
    clusters <- .get_clusters(clus_obj=clus_obj, k=k,
                            per_group=FALSE)
  } else {
    clusters <- clus_obj
  }
  precs <- unlist(precip_list)
  precs <- precs[precs>filt]
  
  lower_lim <- max(filt, 1e-3)
  
  y_seq <- seq(lower_lim, max(precs), length.out=length_out)
  
  if(length(clusters)==length(precip_list)) { #i.e. the partitioning case
    clusters <- rep(clusters,
                    times=sapply(precip_list, \(x) sum(x>filt)))
  }
  
  precs_per_cluster <- split(precs, clusters)
  
  # simple Gamma fit (overall)
  Gamma_all <- .calc_gamma(prec_actual=precs, prec_predict=y_seq,
                           type=type)
  
  # aggregated Gamma fit
  weights <- (table(clusters)/length(clusters)) |> as.vector()#getModel(clus_obj, as.character(k))@size
  Gamma_aggregated <- sapply(precs_per_cluster,
                             .calc_gamma,
                             prec_predict=y_seq,
                             type=type)
  Gamma_aggregated <- colSums(t(Gamma_aggregated)*weights)
  
  data.frame(rr=y_seq,
             simple_Gamma=Gamma_all,
             aggregated_Gamma=Gamma_aggregated)
  
}

#' Helper function to calculate kullback-leibler divergence between an
#' empirical and a theoretical PDF
#' @param empirical numeric vector of empirical PDF
#' @param fitted numeric vector of corresponding theoretical PDF
.get_KL <- function(empirical, fitted) {
  mat <- rbind(empirical, fitted)
  mat <- mat/rowSums(mat) #normalization of each distribution
  philentropy::KL(mat)
}