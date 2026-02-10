#rewriting it to bootstrap events - should be much faster too
#21.07.25

#rewriting flexclust:::MClapply so that it'll handle *actually multiple*
#cores, and not just two (default value of parallel::mclapply)
MClapply <- function (X, FUN, multicore = TRUE, mc.cores=2, ...) 
{
  if (inherits(multicore, "cluster")) 
    parLapply(multicore, X, FUN)
  else if (multicore) 
    mclapply(X, FUN, mc.cores=mc.cores, ...)
  else lapply(X, FUN, ...)
}

#bootstrapping of flexmix at group(=event) level
bootEventmix <- function(x, k, formula=x~1, model=FLXMCnorm1(), filter=0,
                         nboot = 100, correct = TRUE, seed = NULL, multicore = TRUE,
                         mc.cores = 2, verbose = FALSE, ...) {
  MYCALL <- match.call()
  if (!is.null(seed)) 
    set.seed(seed)
  seed <- round(2^31 * runif(nboot, -1, 1))
  nk <- length(k)
  nx <- nrow(x)
  xlong <- unnest(x, precipitation) |> 
    filter(precipitation > filter)
  index1 <- matrix(integer(1), nrow = nx, ncol = nboot)
  index2 <- index1
  for (b in 1:nboot) {
    index1[, b] <- sample(1:nx, nx, replace = TRUE)
    index2[, b] <- sample(1:nx, nx, replace = TRUE)
  }
  BFUN <- function(b) {
    if (verbose) {
      if ((b%%100) == 0) 
        cat("\n")
      if ((b%%10) == 0) 
        cat(b, "")
    }
    set.seed(seed[b])
    
    x1 <- x[index1[, b], , drop = FALSE] |> 
 #     mutate(id=seq_along(id)) |> #done to prevent y~x|id from forcing all 'new' events into one cluster --> BUT causes loads of issues later on. Just leaving it out, see inst/scripts/bootEventmix_outtakes.R
      unnest(precipitation) |> 
      filter(precipitation > filter)
    # events1 <- x$id[index1[,b]]
    
    x2 <- x[index2[, b], , drop = FALSE] |> 
#      mutate(id=seq_along(id)) |>
      unnest(precipitation) |> 
      filter(precipitation > filter)
    # events2 <- x$id[index2[,b]]
    
    s1 <- tryCatch({stepFlexmix(x1, k = k,
                      formula=formula, model=model,
                      verbose = FALSE, drop = FALSE,
                      ...)},
                   error=function(e) {
                     return(e$message)
                   })
    s2 <- tryCatch({stepFlexmix(x2, k = k,
                                formula=formula, model=model,
                                verbose = FALSE, drop = FALSE,
                                ...)},
                   error=function(e) {
                     return(e$message)
                   })
    
    ks <- list(k1=s1, k2=s2)
    dtcls <- sapply(c(s1, s2), data.class)
    
    for(i in 1:2) {
      if(dtcls[[i]]=="stepFlexmix") {
        ks[[i]] <- sapply(ks[[i]]@models, \(y) max(clusters(y)))
      } #ks$k1@k extracts k0 (the input k) not k (the used k)
    }
    
    
    #1) count the frequency of actual k use
    k_counts <- sapply(ks, \(y) {
      sapply(k, \(x) sum(y==x))
    })
    rownames(k_counts) <- k
    tot_k_counts <- rowSums(k_counts) #that's the info I'll extract in the end
    
    if(all(dtcls=="stepFlexmix")) {
      #2) for duplicates, choose the 'better' model
      dups <- which(k_counts>1, arr.ind=TRUE, useNames=TRUE) |> 
        as.data.frame() |> 
        dplyr::mutate(k=k[row], #rownames_to_column no workey, because if a row is duplicate (i.e. both cols, the rownames are destroyed by as.data.frame)
                      best_k0=NA)
      
      for(i in seq_len(nrow(dups))) {
        #find the 'better' model
        h <- which(dups$k[i]==ks[[dups$col[i]]]) |> 
          names()
        if(dups$col[i]==1) md <- s1 else md <- s2
        dups$best_k0[i] <- sapply(h, \(y) ICL(getModel(md, y))) |> 
          which.min() |> names()
        #now filter the ks
        h <- names(which(ks[[dups$col[i]]]==dups$k[i])) |> 
          setdiff(dups$best_k0[i])
        ks[[dups$col[i]]] <- ks[[dups$col[i]]][
          setdiff(names(ks[[dups$col[i]]]), h)
        ]
      }
      
      #3) map the k0_i's by the k values
      ks <- lapply(ks, enframe, name='k0', value='k')
      ks <- dplyr::full_join(ks[[1]], ks[[2]], by='k',
                             suffix=paste0('_', 1:2))
      
      #4) 'stretch' k with NAs #die wieder dazuzutun, ist glaube ich eigentlich unnötig, aber schauen wir später #hi
      ks <- rbind(ks, c(NA, setdiff(k, ks$k), NA)) |> 
        dplyr::arrange(k)
      
      success_pairs <- na.omit(ks)
      
      repeat { #In general, I get that this is some kind of failure of fit check, but what for? Where is it used?
        count <- sapply(seq_len(nrow(success_pairs)), function(i) {
          m1 <- getModel(s1, success_pairs$k0_1[i])
          m2 <- getModel(s2, success_pairs$k0_2[i])
          sapply(list(m1, m2), function(m) {
            any(m@size < 2) || length(m@size) != 
              m@k# || any(m@clusinfo[, 2] == 0) #not sure what the modbased equivalent would be here
          })
        })
        if (!any(count) || 1) 
          break
        rejected <- rejected + 1
        count <- 0
      }
      clust1 <- clust2 <- matrix(integer(1), nrow = nrow(xlong), ncol = nk,
                                 dimnames=list(NULL, k))
      post1 <- post2 <- vector('list', length=nk) |> 
        setNames(k)
      comps1 <- comps2 <- vector('list', length=nk) |> 
        setNames(k)
      rand <- double(nk) |> setNames(k)
      for (l in k) {
        lc <- as.character(l)
        if(!(l %in% success_pairs$k)) {
          clust1[, lc] <- clust2[, lc] <- NA
          post1[[lc]] <- post2[[lc]] <- list(scaled=NA, unscaled=NA)
          comps1[[lc]] <- comps2[[lc]] <- NA
          rand[lc] <- NA
        } else {
          lr <- which(l==success_pairs$k)
          cl1 <- getModel(s1, success_pairs$k0_1[lr])
          cl2 <- getModel(s2, success_pairs$k0_2[lr])
          clust1[, lc] <- clusters(cl1, newdata = as.data.frame(xlong)) #clusters() doesn't like tibbles
          clust2[, lc] <- clusters(cl2, newdata = as.data.frame(xlong))
          post1[[lc]] <- cl1@posterior
          post2[[lc]] <- cl2@posterior
          comps1[[lc]] <- cl1@components
          comps2[[lc]] <- cl2@components
          rand[lc] <- randIndex(table(clust1[, lc], clust2[, lc]),
                                correct = correct)
        }
        #legacy code, what is it for? (I guess I don't need it in my matching setup right?)
        # if (nrow(cl1@centers) < k[l]) {
        #   extra <- matrix(NA, ncol = ncol(cl1@centers), 
        #                   nrow = k[l] - nrow(cl1@centers))
        #   cent1[[l]] <- rbind(cl1@centers, extra)
        # }
        # if (nrow(cl2@centers) < k[l]) {
        #   extra <- matrix(NA, ncol = ncol(cl2@centers), 
        #                   nrow = k[l] - nrow(cl2@centers))
        #   cent2[[l]] <- rbind(cl2@centers, extra)
        # }
      }
      list(post1 = post1, post2 = post2, clust1 = clust1, clust2 = clust2, 
           comps1 = comps1, comps2 = comps2, rand = rand, tot_k_counts = tot_k_counts,
           error=NA)
    } else {
      
      if(any(dtcls)=="stepFlexmix") {
        mds <- c(s1, s2)
        post <- sapply(mds, \(x) {
          tryCatch({sapply(x@models, posterior)}, 
                   error=function(e) NA)
        })
        clust <- sapply(mds, \(x) {
          tryCatch({sapply(x@models, clusters)}, 
                   error=function(e) NA)
        })
        comps <- sapply(mds, \(x) {
          tryCatch({sapply(x@models, \(y) y@components)}, 
                   error=function(e) NA)
        })
        err <- ks[[which(dtcls!="stepFlexmix")]]
      } else {
        post <- clust <- comps <- c(NA, NA)
        err <- paste(ks, collapse="_")
      }
      
      list(post1 = post[[1]], post2 = post[[2]], clust1 = clust[[1]], clust2 = clust[[2]], 
           comps1 = comps[[1]], comps2 = comps[[2]], rand = NA, tot_k_counts = tot_k_counts,
           error=err)
    }
    

  }
  z <- MClapply(as.list(1:nboot), BFUN, multicore = multicore, mc.cores=mc.cores)
  clust1 <- sapply(z, `[[`, 'clust1', simplify='array')
  clust2 <- sapply(z, `[[`, 'clust2', simplify='array')
  
  comps1 <- lapply(z, \(y) sapply(y$comps1, `[[`, 1))
  comps2 <- lapply(z, \(y) sapply(y$comps2, `[[`, 1))
  
  #TODO: praktischere Umformung für @posteriors und @components schreiben.
  
  post1 <- lapply(z, \(y) {
    scaled <- sapply(y$post1, `[[`, 'scaled')
    unscaled <- sapply(y$post1, `[[`, 'unscaled')
    list(scaled=scaled, unscaled=unscaled)
  })
  post2 <- lapply(z, \(y) {
    scaled <- sapply(y$post2, `[[`, 'scaled')
    unscaled <- sapply(y$post2, `[[`, 'unscaled')
    list(scaled=scaled, unscaled=unscaled)
  })
  
  comps1 <- lapply(z, \(y) {
    sapply(y$comps1, `[[`, 1)
  })
  comps2 <- lapply(z, \(y) {
    sapply(y$comps2, `[[`, 1)
  })
  
  error <- sapply(z, `[[`, "error")
  
  if (nk > 1) {
    rand <- t(sapply(z, \(x) x$rand))
    tot_k_counts <- t(sapply(z, \(x) x$tot_k_counts))
  } else {
    rand <- as.matrix(sapply(z, \(x) x$rand))
    tot_k_counts <- as.matrix(sapply(z, \(x) x$tot_k_counts))
  }
  
  
  if (verbose) 
    cat("\n")
  list(k = as.integer(k),
       posteriors1 = post1, posteriors2 = post2,
       cluster1 = clust1, cluster2 = clust2,
       components1 = comps1, components2 = comps2,
       index1 = index1, index2 = index2,
       rand = rand, total_ks_actuallyUsed = tot_k_counts, error = error,
       call = MYCALL)
  #die nette print-Methode der definierten Klasse fehlt halt jetzt
}

if(FALSE) {
  h <- rename(scaled10$`30`, precipitation = rr) |> #not yet unnested
    bootEventmix(k=K, formula=y~x|id, filter=0, #here I'd put 1.27/6 (caution, I'm always filtering 'bigger'. Is 1.27/6 written for bigger equal?)
                 model=FLXMRglm(frm, family='Gamma'),
                 control=list(iter.max=100),
                 nboot=5,
                 multicore=FALSE)
}