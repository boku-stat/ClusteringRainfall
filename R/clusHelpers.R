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
 