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
  
  df = cbind(observed.time, eventind, X)
  return(list(observed.time = observed.time, eventind = eventind))
}
