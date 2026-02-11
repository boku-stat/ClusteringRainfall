#' Wrapper to conduct FAMD analysis and create ggplot colored by
#' clustering vector; and with shape determined by `flash`.
#' @param dat Dataset to be visualized (with only the columns selected for FAMD analysis)
#' @param clusters Vector of partitions (integers; or factor)
#' @param dims Dimensions to be used for visualization. List of vectors for x and y
#' @param scale_factor Factor by which the visual representation of the original variables will be scaled
#' @param title optional plot title
#' @param subtitle optional plot subtitle
#' @param varlab_length number of characters to be used for variable labels.
#'        If `NULL`, the full variable names are used.
#' @param varlab_abbrev_show show the legend to the abbreviated variables? `TRUE/FALSE`.
#'        Default: `TRUE`
#' @param varlab_size variable label size; default 3.5
#' @param varlab_nudge factor by which variable label is nudged; default 1.1
#' @param varlab_color char/hexcode with desired color for the variable labels. Default: "blue"
#' @param vararrow_color char/hexcode with desired color for the variable arrow. Default: "blue"
plt_famd <- function(dat, clusters, dims=list(x=2:4, y=1),
                     scale_factor=5, title=NULL, subtitle=NULL,
                     varlab_length=NULL, varlab_abbrev_show=TRUE,
                     varlab_size=3.5, varlab_nudge=1.1,
                     varlab_color='blue', vararrow_color='blue'){
  # Factorize flash if that hasn't been done yet
  if(!is.factor(dat$flash)) dat$flash <- factor(dat$flash)
  dim <- lapply(dims, \(d) sprintf("Dim.%d", d))
  famd_result <- FactoMineR::FAMD(dat, graph = FALSE)
  
  # Extract coordinates for all dimensions, and attach cluster and flash info
  ind_coords <- famd_result$ind$coord |>
    as_tibble() |> 
    mutate(cluster = factor(clusters),
           flash = dat$flash) |> 
    pivot_longer(all_of(dim$x), names_to = "x",
                 values_to = "values_x") |> 
    pivot_longer(all_of(dim$y), names_to = "y",
                 values_to = "values_y")
  # Extract variable coordinates and scale them
  var_coords <- famd_result$var$coord |>
    as_tibble() |> 
    mutate(variable = rownames(famd_result$var$coord),
           varlab = variable) |> 
    mutate(across(where(is.numeric), ~ .x * scale_factor)) |> 
    pivot_longer(all_of(dim$x), names_to = "x",
                 values_to = "values_x") |> 
    pivot_longer(all_of(dim$y), names_to = "y",
                 values_to = "values_y")
 
  abbrev_caption <- NULL
  if(!is.null(varlab_length)) { #this is not yet generalized. Confusion would arise if several variables had the same first letter(s)
    if(!is.numeric(varlab_length) || length(varlab_length) != 1 || varlab_length < 1)
      stop("varlab_length must be a single numeric value >= 1")
    var_coords$varlab <- substr(var_coords$varlab, 1, varlab_length) |> 
      str_to_title()
    if(varlab_abbrev_show) {
      abbrev_caption <- select(var_coords, varlab, variable) |> 
        unique() |> sort() |> 
        apply(1, paste, collapse=' ... ') |> 
        paste(collapse='; ')
    }
  }
  # Extract variance percentage
  labs <- lapply(dims, \(d) {
    famd_result$eig[,"percentage of variance"][d] |> 
      round(2)
  })
  labs <- sapply(names(dim), \(d) {
    sprintf("%s (%g%%)", dim[[d]], labs[[d]]) |> 
      setNames(dim[[d]])
  }, simplify = FALSE)

    ggplot() +
      geom_point(data = ind_coords,
                 aes(x=values_x, y=values_y,
                     color=cluster, shape=flash), size=2) +
      geom_segment(data = var_coords, x=0, y=0,
                   aes(xend = values_x, yend = values_y),
                   arrow = arrow(length = unit(0.2, "cm")), color = vararrow_color) +
      geom_text(data = var_coords, 
                aes(x = values_x*varlab_nudge, y = values_y*varlab_nudge,
                    label = varlab),
                color = varlab_color, size=varlab_size) +
      labs(x=NULL, y=NULL, title=title, subtitle=subtitle,
           caption=abbrev_caption) +
      facet_grid(y~x, scales="free",
                 labeller=labeller(x=labs$x, y=labs$y),
                 switch="both") +
      theme_minimal() +
      theme(strip.placement = "outside")
  
}

if(FALSE) {
  #Usage
  plt_famd(dat=x1, clusters=pam_result_k3$clustering)
  plt_famd(dat=x1, clusters=flx.test@cluster)
  plt_famd(dat=x1, clusters=flx.gauss@cluster)
  plt_famd(dat=x1, clusters=m@cluster)
}
#Source: adapted from github.com/glaaha/mixrain/R/famd_analysis.R



#' Plot Partial Dependency Plots of a `flexmix` clusterwise regression object: Step 1: obtain data for plotting
#' @param model object of class `flexmix` that has been calculated with covariates (i.e. clusterwise regression)
#' @param effect name of the covariate whose partial dependency is to be plotted (character variable, length=1)
#' @param mean_vars vector of the names of the variables whose influence shall be averaged (character variable, length>=1)
#' @param fac_var name of the factor covariate whose influence shall also be visualized (character variable, length=1)
#'                Note that this function is only adapted to work with factors with `contr.treatment` type contrasts.
.get_pdp_data <- function(model, effect, mean_vars, fac_var)
{
  x <- model@model[[1]]@x 
  out <- as_tibble(x)
  cl <- model@cluster
  out$cluster <- cl
  cnames <- colnames(x)
  #reconstruct fac_var
  fclvs <- paste0('^', fac_var) |> 
    grep(cnames, value=TRUE)
  lvs <- gsub(fac_var, '', fclvs)
  fac_var <- out[fclvs]
  if(any(grepl('Intercept', cnames))) {
    lvs <- c(lvs, 'reference_group')
    fac_var$reference_group <- as.numeric(!rowSums(fac_var))
  }
  out$fac_var <- fac_var <- apply(fac_var, 1, \(x) lvs[as.logical(t(x))])

  ind <- grep(paste(mean_vars, collapse="|"), x=cnames)
  # Compute means
  out <- out |> group_by(cluster, fac_var) |> 
    mutate(across(ind, mean)) |> ungroup()
  # shapes <- lapply(model@components, function(x) x[[1]]@parameters$shape) #what are these values created for? They're not used again
  # cf <- lapply(model@components, function(x) x[[1]]@parameters$coef)
  k <- model@k
  # Add column for prediction
  out$predictions <- NA
  for(i in 1:k)
  {
    xx <- out[cl==i,cnames] |> as.matrix()
    fpred <- model@components[[i]][[1]]@predict
    pred <- fpred(xx)[,1]
    out$predictions[out$cluster == i] <- pred
  }
  out <- out[,grep('Intercept', colnames(out), invert=TRUE)] # Do not need Intercept column for final output #caution, in our preferred model specification we've turned off the Intercept
  out
}
#Source: adapted from N: Projekte/EROSA-Stat/Tawes_stations/pdp_flexmix.R

#' Plot Partial Dependency Plots of a `flexmix` clusterwise regression object: Step 2: plot
#' @param model object of class `flexmix` that has been calculated with covariates (i.e. clusterwise regression)
#' @param eval_vars vector of the names of the variables which shall be evaluated, as effect and/or mean variable
#' @param effect_vars vector of the names of the variables whose effect shall be investigated. Default: equal to `eval_vars`.
#' @param effect_type_main Logical. Shall the main effect of `effect_vars` be investigated? (Default=TRUE. FALSE --> interaction effects are investigated)
#' @param fac_var name of the factor covariate whose influence shall also be visualized (character variable, length=1)
#' @param title optional plot title
#' @param subtitle optional plot subtitle
#' Note that this function is only adapted to work with factors with `contr.treatment` type contrasts.
plt_dpd_flexmix <- function(model, eval_vars,
                            effect_vars=eval_vars,
                            effect_type_main=TRUE,
                            fac_var,
                            title=NULL, subtitle=NULL) {

  df <- sapply(effect_vars, \(x) {
    mean_vars <- setdiff(eval_vars, x)
    df <- .get_pdp_data(model=model, effect=x,
                        mean_vars=mean_vars, fac_var=fac_var)
    effect_col <- grep(x, colnames(df), value=TRUE) |> 
      grep(':', x=_, value=TRUE, invert=effect_type_main)
    #if(grepl('log', effect_col)) #I *think* I need to do this always, cuz after all, it's a Gamma GLM with a log link
    df$value <- exp(df[[effect_col]])
    df$name <- effect_col
    df
  }, simplify = FALSE) |> 
    bind_rows(.id="effect")
  
  lbl <- function(x) {
    gsub("_", " ", x) |> 
      str_to_title()
  }
  
  ggplot(df, aes(x=value, y=predictions - 10,
                 col=as.factor(cluster),
                 linetype=as.factor(fac_var))) +
    geom_line(linewidth=1.2) +
    facet_wrap(~effect, scales="free",
               strip.position="bottom",
               labeller=as_labeller(lbl)) +
    labs(x=NULL, y="Precipitation (mm)",
         colour="Cluster", linetype=str_to_title(fac_var),
         title=title, subtitle=subtitle) +
    theme_bw() + 
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          panel.grid = element_line(linewidth = .25), 
          legend.position = "bottom")
  
}



#' Helper that merges cluster vector to data and then unnests and filters
#' @param dat Dataset containing nested precipitation column called 'rr'
#' @param clusters vector of event-level clusters
#' @param filt filter value for the precipitation values (default=0)
.mod_dat <- function(dat, clusters, filt) {
  mutate(dat, clusters=as.factor(clusters)) |> 
    unnest(rr) |> 
    filter(rr > filt)
}

#' Wrapper for ECDF plotting - both by group (=cluster) and total
#' @param dat Dataset to be visualized - at event level, contains nested precipitation column (called "rr")
#' @param clus_obj One of:
#'                  - a clustering object of class `stepFlexclust` or `stepFlexmix`,
#'                  - vector of cluster assignments,
#'                  - list of cluster assignment vectors, named.
#' @param k the desired number of clusters out of the `k` values tested in the clustering object
#'          (only used if `clus_obj` is of class `stepFlexmix`/`stepFlexclust`)
#' @param filt filter value for the precipitation values (default=0)
#' @param title plot title (character vector, length=1)
#' @param subtitle Character vector of plot subtitle(s). If length>1, the vector must
#'                 named the same as the pertaining `clus_obj` list items.
plt_ecdf <- function(dat, clus_obj, k=NULL, filt=0,
                          title=NULL, subtitle=NULL) {

  if(data.class(clus_obj) %in% c('stepFlexclust', 'stepFlexmix')) {
    if(is.null(k)) stop('Please specify desired number of clusters k!')
    clusters <- .get_clusters(clus_obj=clus_obj, k=k,
                              per_group=FALSE)
    df <- .mod_dat(dat=dat, clusters=clusters, filt=filt) |> 
      mutate(model=data.class(clus_obj))
    labs <- list(model=subtitle)
  } else if(data.class(clus_obj)=='list') {
    k <- max(unlist(clus_obj)) #only for the number of colors in the palette
    df <- lapply(clus_obj,
                 \(cl) .mod_dat(dat=dat, clusters=cl, filt=filt)) |> 
      bind_rows(.id='model')
    labs <- subtitle
  } else {
    k <- max(clus_obj)
    df <- .mod_dat(da=dat, clusters=clus_obj, filt=filt) |> 
      mutate(model='model')
    labs <- list(model=subtitle)
  }
  
  ggplot(df, aes(x=rr)) +
    stat_ecdf(geom='step', aes(color='all observations')) +
    stat_ecdf(geom='step', aes(color=clusters)) +
    facet_wrap(~model, scales='free', dir="rt", #quick dirty hack to get 'P' before 'C'
               labeller=as_labeller(labs)) +
    scale_x_log10() +
    theme_bw() +
    labs(title=title,
         x='Precipitation (mm)', y='ECDF') +
    scale_color_manual(name='Group',
                       values=c(`all observations` = 'black',
                                setNames(scales::hue_pal()(k),
                                         1:k)),
                       labels = \(x) ifelse(x == 'all observations', x,
                                            paste('Cluster', x))) +
    theme(strip.background = element_blank(),
          legend.position = 'bottom')
}
#Source: adapted from github.com/glaaha/mixrain/R/dist_analysis.R

#' Visual comparison of overall Gamma vs. clusterwise-aggregated Gamma distributions
#' to precipitation ECDF
#' @param dat Dataset to be visualized that contains rr (=precipitation) column nested at event level
#' @param clus_vec One of:
#'                  - event level vector of clusters,
#'                  - list of event level vectors of clusters.
#' @param filt lower limit for (E)CDF fitting of precipitation values
#' @param title plot title (character vector, length=1)
#' @param subtitle Character vector of plot subtitle(s). If length>1, the vector must
#'                 named the same as the pertaining `clus_obj` list items.
plt_cdf_comp <- function(dat, clus_vec, filt, title=NULL, subtitle=NULL) {

  lvs <- data.frame(
    levels=c('Empirical', 'simple_Gamma', 'aggregated_Gamma'),
    labels=c('Empirical CDF',
             'Gamma (overall fit)',
             'Gamma (cluster-aggregated)'),
    colors=c('black', 'grey', 'firebrick1')
  )
  
  if(data.class(clus_vec)=='list') {
    ecdf <- lapply(clus_vec, \(x) {
      aggregate_dist_calc(dat$rr, x,
                          filt=filt, type='ECDF') |> 
        pivot_longer(matches('Gamma')) |> 
        mutate(name=factor(name,
                           levels=lvs$levels))
    }) |> bind_rows(.id='model')
    labs <- subtitle
  } else {
    ecdf <- aggregate_dist_calc(dat$rr, clus_vec,
                                filt=filt, type='ECDF') |> 
      pivot_longer(matches('Gamma')) |> 
      mutate(name=factor(name,
                         levels=lvs$levels),
             model='model')
    labs <- list(model=subtitle)
  }
  
  ggplot() +
    stat_ecdf(data=filter(unnest(dat, rr), rr>filt), 
              aes(x=rr, col='Empirical'),
              geom='step') +
    geom_line(data=ecdf,
              aes(x=rr, y=value, col=name)) +
    facet_wrap(~model, scales='free', dir="rt", #quick dirty hack to get 'P' before 'C'
               labeller=as_labeller(labs)) +
    scale_color_manual(
      name = 'Distribution',
      values = deframe(select(lvs, -labels)),
      labels = deframe(select(lvs, -colors)),
      breaks = lvs$levels
    ) +
    scale_x_log10(limits=c(0.05, NA)) +
    labs(title=title,
         x='Precipitation (mm)', y='Cumulative Distribution Function') +
    theme_bw() +
    theme(strip.background = element_blank(),
          legend.position = 'bottom')
}