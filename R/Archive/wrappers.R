library(cmdstanr)
library(posterior)
library(bayesplot)

## directory where stan files are located
standir <- '/proj/ibrahimlab/strapp_paper3/Stan'

hist_post_stan            <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_hist_post.stan'), include_paths = standir)
pwe_hist_post_stan        <- cmdstanr::cmdstan_model(file.path(standir, 'pwe_logistic_hist_post.stan'), include_paths = standir)

curepwe_logistic_mix_stan           <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_mixture.stan'), include_paths = standir)
curepwe_logistic_mix_strapp_stan    <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_mixture_strapp.stan'), include_paths = standir)
curepwe_logistic_mix_genstrapp_stan <- cmdstanr::cmdstan_model(file.path(standir, 'curepwe_logistic_mixture_genstrapp.stan'), include_paths = standir)


#' Sample from the posterior of the logistic historical data (current data model: curePWE)
#' 
#' @param fmla.hist formula for historical data
#' @param histdata a data frame for historical (glm) data
#' @param a0
#' @param pc_upper upper bound of the uniform prior on p_cured if isBeta = 0
#' @param isBeta indicator of whether a Beta prior on p_cured should be used (1 if using a Beta prior)
#' @param pc_shape the shape/rate parameter of the Beta prior, i.e., p_cured ~ Beta(pc_shape, pc_shape)
#' @return object of type cmdstanrfit

hist_post <- function(
  fmla.hist, 
  histdata,
  a0 = 0.2, 
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4,
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
    , 'pc_upper' = pc_upper
    , 'isBeta' = isBeta
    , 'pc_shape' = pc_shape
  )
  # return(standat)
  smpl <- hist_post_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_hist_post"
  smpl
}


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
  # return(standat)
  smpl <- pwe_hist_post_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "pwe_logistic_hist_post"
  smpl
}


#' Sample from curepwe_logistic_mixture
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

curepwe_mix <- function(
  fmla.cur, data,
  breaks,
  pc_upper = 1,
  isBeta = 1,
  pc_shape = 4,
  weights,
  means,
  covars
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
  )
  # return(standat)
  smpl <- curepwe_logistic_mix_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_mixture"
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
  )
  # return(standat)
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
  )
  # return(standat)
  smpl <- curepwe_logistic_mix_genstrapp_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "curepwe_logistic_mixture_genstrapp"
  smpl
}


