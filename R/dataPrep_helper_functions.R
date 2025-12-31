#### Helper functions

get_estimate <- function(x)
{
  m <- mledist(data = x, distr = "gamma")
  m$estimate
}

get_prob_cluster <- function(x, model, var_name = "rr")
{
  x <- x %>% mutate(cl = model@cluster)
  x <- x[c("cl", var_name)] %>% rename(y := !!as.name(var_name))
  ncl <- unique(x$cl)
  par <- lapply(ncl, function(cc) get_estimate(x %>% filter(cl == cc) %>% pull(y)))
  # Get estimate for full time series
  par_full <- get_estimate(x$y)
  x <- x %>% mutate(fp = pgamma(y, shape = par_full[1], rate = par_full[2]))
  # Get estimate for empirical distribution
  x <- x %>% mutate(ep = rank(y)/(n() + 1)) # Weibull plotting position 
  # Get estimate for each cluster
  x <- x %>% mutate(clp = NA)
  for(i in 1:length(ncl))
  {
    ii <- which(x$cl == ncl[i])
    x$clp[ii] <- pgamma(x$y[ii], shape = par[[i]][1], rate = par[[i]][2])
  }
  # Get weights
  n <- length(x$cl)
  w <- lapply(ncl, function(cc) length(x$cl[x$cl == cc])/n)
  l <- length(ncl)
  cp <- sapply(seq(l), function(ind) w[[ind]]*pgamma(x$y, shape = par[[ind]][1], rate = par[[ind]][2]))
  x <- x %>% mutate(ap = apply(cp, 1, sum))
  x
}

# This function would add weighted samples!

get_sample <- function(x)
{
  x <- x[x>0]
  if(length(x) == 1)
  {
    out <- x
  }
  else
  {
    out <- sample(x, size = 1, prob = x/sum(x))
  }
  out
}

# Get precipitation data
get_precip <- function(start, end, df)
{
  rr <- df %>% filter(between(time, start, end)) %>% pull(rr)
  rr
}

# Get flash data
get_flash_data <- function(start, end, flash)
{
  ff <- flash %>% filter(between(date, start - hours(1), end)) %>% pull(flash)
  res <- ifelse(sum(ff) > 0, 1, 0)
  res
}

get_de <- function(x)
{
  epdf <- diff(c(0,sort(x)))
  epdf/sum(epdf)
}

# get_KL <- function(x,y) # it is assumed that x and y are both cumulative distribution functions!
# {
#   xdf <- get_de(x = x)
#   ydf <- get_de(x = y)
#   m <- rbind(xdf, ydf)
#   philentropy::KL(m)
# }
# 
# back_transform <- function(x, y)
# {
#   maxy <- max(y)
#   miny <- min(y)
#   xtrans <- (x - 0.01)*(maxy - miny) + miny
#   xtrans
# }
