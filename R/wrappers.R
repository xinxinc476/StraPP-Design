library(cmdstanr)
library(posterior)
library(bayesplot)

## directory where stan files are located
standir <- '/proj/ibrahimlab/strapp_paper3/Stan'

pwe_hist_post_stan        <- cmdstanr::cmdstan_model(file.path(standir, 'pwe_logistic_hist_post.stan'), include_paths = standir)

curepwe_logistic_mix_strapp_stan    <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_mixture_strapp.stan'), include_paths = standir)
curepwe_logistic_mix_genstrapp_stan <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_mixture_genstrapp.stan'), include_paths = standir)
curepwe_logistic_pp_stan            <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_pp.stan'), include_paths = standir)
curepwe_refprior_stan               <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_refprior.stan'), include_paths = standir)

#' Sample from the posterior of the logistic historical data (current data model: PWE)
#' 
#' @param fmla.hist formula for historical data
#' @param histdata a data frame for historical (glm) data
#' @param a0
#' @return object of type cmdstanrfit

pwe_hist_post <- function(
  fmla.hist, 
  histdata,
  a0 = 0.2, 
  ...
) {
  
  y0name <- all.vars(fmla.hist)[1]
  y0     <- histdata[, y0name]
  X0     <- model.matrix(fmla.hist, histdata)
  
  ## Make sure no design matrices have intercepts
  if ( '(Intercept)' %in% colnames(X0) )
    X0 <- X0[, -1, drop = FALSE]
  
  standat <- list(
    'n0'    = nrow(X0)
    , 'p'   = ncol(X0)
    , 'y0'  = y0
    , 'X0'  = X0
    , 'a0'  = a0
  )
  
  smpl <- pwe_hist_post_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "pwe_logistic_hist_post"
  smpl
}


#' Sample from curepwe_logistic_mixture_strapp
#' new approximation
#' 
#' @param fmla.cur formula for current data
#' @param data a data frame for current (survival) data
#' @param breaks
#' @param weights a vector of length K, mixing proportions for the Gaussian mixture prior
#' @param means a p by K mean matrix, where p is the number of covariates, and K is the number of components
#' means[, i] is the mean vector of the ith component in the Gaussian mixture model 
#' @param covars a p by p by K covariance array
#' covars[, , i] is the covariance matrix of the ith component in the Gaussian mixture model 
#' @return object of type cmdstanrfit

curepwe_mix_strapp <- function(
  fmla.cur, data,
  breaks,
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4,
  weights,
  means,
  covars
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
  , ...
) {
  
  ## Extract names and variables for response, censoring, etc.
  y1name    <- all.vars(fmla.cur)[1]
  eventname <- all.vars(fmla.cur)[2]
  
  y1       <- data[, y1name]
  eventind <- data[, eventname]
  
  X1 <- model.matrix(fmla.cur, data)
  if ( '(Intercept)' %in% colnames(X1) )
    X1 <- X1[, -1, drop = FALSE]
  
  p  <- ncol(X1)
  
  ## Shape and rate for hazard as vector
  J <- length(breaks) - 1  ## number of intervals
  
  ## Create index giving interval into which obs failed / was censored
  intindx <- rep(NA, nrow(X1))
  for ( j in 1:J ) {
    intindx[ y1 >= breaks[j] & ( y1 <= breaks[j+1] )  ] <- j
  }
  
  K <- length(weights)
  
  # ensure that the mean matrix has dimension = c(p, K)
  means <- matrix(means, nrow = p, ncol = K)
  
  # reshape the covariance array to have dimension = c(K, p, p)
  if ( (length(covars) == 1) & (K > 1) ) {
    covars <- rep(covars, K)
  }
  if ( is.null(dim(covars)) ){
    covars <- array(covars, dim = c(1, 1, K))
  }
  if ( is.na(dim(covars)[3]) ){
    covars <- array(covars, dim = c(p, p, K))
  }
  covars_reshaped <- array(dim = c(K, p, p))
  for(i in 1:K) {
    covars_reshaped[i, , ] = covars[ , , i]
  }
  
  standat <- list(
      'n1'           = nrow(X1)
    , 'J'            = J
    , 'p'            = p
    , 'y1'           = y1
    , 'X1'           = X1
    , 'intindx'      = intindx
    , 'death_ind'    = eventind
    , 'breaks'       = breaks
    , 'pc_upper'     = pc_upper
    , 'isBeta'       = isBeta
    , 'pc_shape'     = pc_shape
    , 'K'            = K
    , 'weights'      = weights
    , 'means'        = means
    , 'covars'       = covars_reshaped
    , 'rel_tol'      = rel.tol
    , 'f_tol'        = f.tol
    , 'max_steps'    = max.steps
  )
  
  smpl <- curepwe_logistic_mix_strapp_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_mixture_strapp"
  smpl
}


#' Sample from curepwe_logistic_mixture_genstrapp
#' new approximation
#' 
#' @param fmla.hist formula for historical data
#' @param fmla.cur formula for current data
#' @param histdata data frame for historical (logistic) data
#' @param data a data frame for current (survival) data
#' @param breaks
#' @param weights a vector of length K, mixing proportions for the Gaussian mixture prior
#' @param means a p by K mean matrix, where p is the number of covariates, and K is the number of components
#' means[, i] is the mean vector of the ith component in the Gaussian mixture model 
#' @param covars a p by p by K covariance array
#' covars[, , i] is the covariance matrix of the ith component in the Gaussian mixture model 
#' @return object of type cmdstanrfit

curepwe_mix_genstrapp <- function(
  fmla.hist, fmla.cur, 
  histdata, data,
  breaks,
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4,
  weights,
  means,
  covars
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
  , ...
) {
  
  ## Extract names and variables for response, censoring, etc.
  y1name    <- all.vars(fmla.cur)[1]
  eventname <- all.vars(fmla.cur)[2]
  
  y1       <- data[, y1name]
  eventind <- data[, eventname]
  
  X1 <- model.matrix(fmla.cur, data)
  if ( '(Intercept)' %in% colnames(X1) )
    X1 <- X1[, -1, drop = FALSE]
  
  X0 <- model.matrix(fmla.hist, histdata)
  if ( '(Intercept)' %in% colnames(X0) )
    X0 <- X0[, -1, drop = FALSE]
  
  p  <- ncol(X1)
  
  ## Shape and rate for hazard as vector
  J <- length(breaks) - 1  ## number of intervals
  
  ## Create index giving interval into which obs failed / was censored
  intindx <- rep(NA, nrow(X1))
  for ( j in 1:J ) {
    intindx[ y1 >= breaks[j] & ( y1 <= breaks[j+1] )  ] <- j
  }
  
  K <- length(weights)
  
  # ensure that the mean matrix has dimension = c(p, K)
  means <- matrix(means, nrow = p, ncol = K)
  
  # reshape the covariance array to have dimension = c(K, p, p)
  if ( (length(covars) == 1) & (K > 1) ) {
    covars <- rep(covars, K)
  }
  if ( is.null(dim(covars)) ){
    covars <- array(covars, dim = c(1, 1, K))
  }
  if ( is.na(dim(covars)[3]) ){
    covars <- array(covars, dim = c(p, p, K))
  }
  covars_reshaped <- array(dim = c(K, p, p))
  for(i in 1:K) {
    covars_reshaped[i, , ] = covars[ , , i]
  }
  
  standat <- list(
      'n1'           = nrow(X1)
    , 'n0'           = nrow(X0)
    , 'J'            = J
    , 'p'            = p
    , 'y1'           = y1
    , 'X1'           = X1
    , 'X0'           = X0
    , 'intindx'      = intindx
    , 'death_ind'    = eventind
    , 'breaks'       = breaks
    , 'pc_upper'     = pc_upper
    , 'isBeta'       = isBeta
    , 'pc_shape'     = pc_shape
    , 'K'            = K
    , 'weights'      = weights
    , 'means'        = means
    , 'covars'       = covars_reshaped
    , 'rel_tol'         = rel.tol
    , 'f_tol'           = f.tol
    , 'max_steps'       = max.steps
  )
  
  smpl <- curepwe_logistic_mix_genstrapp_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_mixture_genstrapp"
  smpl
}


#' Function to sample from a mixture cure rate model where the density for the non-cured population follows a exponential baseline hazard model
#'
#' @param pcure probability of cure
#' @param X
#' @param beta 
#' @param lambda baseline hazard, solved based on 28% survival rates at 5 years (from the 2nd reference in E1684 paper) 
#' @param accrual.max maximum accrual time 
#' For E1690 study, the accrual time was from February 1991 to June 1, 1995 (52/12 years),
#' and the presently analyzed database was current to September 1998
#' @param min.follow.up minimum follow up time for each subject, default is 3
#' @param censor.rate censoring rate, default is log(0.95)/(-5)
#' @return list(observed time, indicator for whether failure observed before censoring) 
rcureexp <- function(pcure, X, beta,
                     lambda = -log(0.28)/5, # solved from exp{-5x}=0.28
                     accrual.max = 52/12,
                     min.follow.up = 3,
                     censor.rate = -log(0.95)/5 # solved from exp{-5x}=0.95, 5% subjects be censored at 5 years
){
  X        <- as.matrix(X)
  n        <- nrow(X)
  eta      <- X %*% beta
  # generate accrual time for each subject
  accrual.time <- runif(n, min = 0, max = accrual.max)
  # administered censor time (the last subject enrolled in study will be observed for min.follow.up years)
  elapsed.time.max <- max(accrual.time) + min.follow.up
  # generate censoring time using exponential distribution with rate = censor.rate
  censor.time <- rexp(n, rate = censor.rate)
  # generate an indicator vector for whether the ith subject is cured
  curedind <- rbinom(n, 1, pcure)
  
  ##
  ## generate failure/event time t from a exponential baseline hazard model
  ##
  t <- rexp( n, rate = exp( log(lambda) + eta ) )
  t[curedind == 1] <- Inf
  
  ## complete data observation time and event indicator
  observed.time <- pmin(t, censor.time)
  eventind      <- as.numeric(t <= censor.time)
  
  ## the complete data elapsed time
  elapsed.time <- accrual.time + observed.time
  
  ## if the elapsed time for a subject is larger than the administered censor time,
  ## then the observed time is set to be (administered censor time - accrual time) and the event indicator is set to 0 
  idx = (elapsed.time > elapsed.time.max)
  observed.time[idx] = elapsed.time.max - (accrual.time[idx])
  eventind[idx] = 0
  
  df           = cbind(observed.time, eventind, X)
  if(ncol(df) == 3) {
    colnames(df) = c('failtime', 'rfscens', 'trt') # treatment only model
  } else{
    colnames(df) = c('failtime', 'rfscens', 'trt', 'age_bin', 'sex', 'node_bin')
  }
  df           = as.data.frame(df)
  return(df)
}


#' Function to sample from curepwe_logistic_pp
curepwe_logistic_pp <- function(
  fmla.cur, fmla.hist, data, histdata,
  breaks, a0 = 0.25, 
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
  , ...
) {
  
  ## Extract names and variables for response, censoring, etc.
  y1name    <- all.vars(fmla.cur)[1]
  eventname <- all.vars(fmla.cur)[2]
  y1        <- data[, y1name]
  eventind  <- data[, eventname]
  X1        <- model.matrix(fmla.cur, data)
  
  if ( is.null(histdata) ) {
    y0name <- ''
    y0 <- 0
    X0 <- matrix(0, 2, 2)
  } else {
    y0name <- all.vars(fmla.hist)[1]
    X0     <- model.matrix(fmla.hist, histdata)
    y0     <- histdata[, y0name]
  }
  
  if ( '(Intercept)' %in% colnames(X1) )
    X1 <- X1[, -1]
  if ( '(Intercept)' %in% colnames(X0) )
    X0 <- X0[, -1]
  
  J <- length(breaks) - 1  ## number of intervals
  
  ## Create index giving interval into which obs failed / was censored
  intindx <- rep(NA, nrow(X1))
  for ( j in 1:J ) {
    intindx[ y1 >= breaks[j] & ( y1 <= breaks[j+1] )  ] <- j
  }
  
  standat <- list(
      'n1'           = nrow(X1)
    , 'n0'           = nrow(X0)
    , 'J'            = J
    , 'p'            = ncol(X0)
    , 'y1'           = y1
    , 'y0'           = y0
    , 'X1'           = X1
    , 'X0'           = X0
    , 'intindx'      = intindx
    , 'death_ind'    = eventind
    , 'breaks'       = breaks
    , 'a0'           = a0
    , 'pc_upper'     = pc_upper
    , 'isBeta'       = isBeta
    , 'pc_shape'     = pc_shape
    , 'rel_tol'      = rel.tol
    , 'f_tol'        = f.tol
    , 'max_steps'    = max.steps
  )
  
  smpl <- curepwe_logistic_pp_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_pp"
  smpl
}


#' Function to sample from curepwe_refprior
curepwe_refprior <- function(
  fmla.cur, data, 
  breaks,  
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4
  , rel.tol = 1e-6, f.tol = 1e-6, max.steps = 1000
  , ...
) {
  
  ## Extract names and variables for response, censoring, etc.
  y1name    <- all.vars(fmla.cur)[1]
  eventname <- all.vars(fmla.cur)[2]
  y1        <- data[, y1name]
  eventind  <- data[, eventname]
  X1        <- model.matrix(fmla.cur, data)
  
  if ( '(Intercept)' %in% colnames(X1) )
    X1 <- X1[, -1, drop = F]
  
  J <- length(breaks) - 1  ## number of intervals
  
  ## Create index giving interval into which obs failed / was censored
  intindx <- rep(NA, nrow(X1))
  for ( j in 1:J ) {
    intindx[ y1 >= breaks[j] & ( y1 <= breaks[j+1] )  ] <- j
  }
  
  standat <- list(
      'n1'           = nrow(X1)
    , 'J'            = J
    , 'p'            = ncol(X1)
    , 'y1'           = y1
    , 'X1'           = X1
    , 'intindx'      = intindx
    , 'death_ind'    = eventind
    , 'breaks'       = breaks
    , 'pc_upper'     = pc_upper
    , 'isBeta'       = isBeta
    , 'pc_shape'     = pc_shape
    , 'rel_tol'      = rel.tol
    , 'f_tol'        = f.tol
    , 'max_steps'    = max.steps
  )
  
  smpl <- curepwe_refprior_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_refprior"
  smpl
}

