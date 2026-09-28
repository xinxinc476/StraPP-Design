#' Inverse of cumulative hazard
#' 
#' @param x vector of cumulative hazards
#' @param lambda J-dimensional vector of hazards
#' @param endpoints (J+1)-dimensional vector of endpoints
#' 
#' @return vector of failure times
inv_cumhaz_pwe <- function(x, lambda, endpoints) {
  ## Get number of hazards / intervals
  J = length(lambda)
  ## Compute cumulative hazard for each interval
  cumhaz = lambda * cumsum(endpoints[2:(J+1)] - endpoints[1:J])
  cumhaz = c(0, cumhaz)
  ## Initialize result
  res <- numeric(length(x))
  
  ## Loop through intervals and compute failure time based on each interval
  for ( j in 1:J ) {
    indx_j      <- x >= cumhaz[j] & x <= cumhaz[j+1]  ## extract obs who fail in specific interval
    if ( length(indx_j) == 0 )
      next
    res[indx_j] <- endpoints[j] + (x[indx_j] - cumhaz[j]) / lambda[j] 
  }
  ## Return result
  res
}


#' Cumulative hazard for a piecewise exponential model
#' 
#' @param t a vector of failure times
#' @param lambda a J-dimensional vector of hazards for each interval
#' @param endpoints a (J+1)-dimensional vector of endpoints
#' 
#' @return vector of cumulative hazards evaluated at each `t`
cumhaz_pwe <- function(t, lambda, endpoints) {
  ## Get number of intervals
  J = length(lambda)
  ## Compute cumulative hazard for each interval
  cumhaz = lambda * cumsum(endpoints[2:(J+1)] - endpoints[1:J])
  cumhaz = c(0, cumhaz)
  ## Initialize result
  res <- numeric(length(t))
  ## Loop through times in each interval and compute cumulative hazard
  for ( j in 1:J ) {
    indx_j      <- (t > endpoints[j]) & (t <= endpoints[j+1])
    if ( length(indx_j) == 0 )
      next
    res[indx_j] <- cumhaz[j] + lambda[j] * (t[indx_j] - endpoints[j])
  }
  ## Return result
  res
}


#' Function to sample from a mixture cure rate model where the density for the non-cured population follows a PWE PH model
#'
#' @param pcure probability of cure
#' @param X
#' @param beta 
#' @param lambda baseline hazard (use J = 8 components)
#' @param endpoints a J+1 vector of endpoints for curePWE intervals, s_1 = 0, s_(J+1) = infinity
#' lambda and endpoints are estimated from fitting PWE model on real dataset (E1690) in SAS
#' @param accrual.max maximum accrual time 
#' For E1690 study, the accrual time was from February 1991 to June 1, 1995 (52/12 years),
#' and the presently analyzed database was current to September 1998
#' @param min.follow.up minimum follow up time for each subject
#' @param censor.rate censoring rate
#'  
#' @return list(observed time, indicator for whether failure observed before censoring) 
rcurepwe <- function(pcure, X, beta,
                     lambda,
                     endpoints,
                     accrual.max = 52/12,
                     min.follow.up = 3,
                     censor.rate = log(0.95)/(-5) # solved from exp{-5x}=0.95, 5% subjects be censored by the end of 5 years
                     ){
  n        <- nrow(as.matrix(X))
  eta      <- X %*% beta
  # generate accrual time for each subject
  accrual.time <- runif(n, min = 0, max = accrual.max)
  # administered censor time (the last subject enrolled in study will be observed for min.follow.up years)
  elapsed.time.max <- max(accrual.time) + min.follow.up
  # generate censoring time using exponential distribution with rate = censor.rate
  censor.time <- rexp(n, rate = censor.rate)
  # generate an indicator vector for whether the ith subject is cured
  curedind <- rbinom(n, 1, pcure)
  
  # generate failure/event time t from a PWE PH model
  z                <- rexp(n, rate = 1)
  x                <- exp(log(z) - eta) # cumulative hazards
  t                <- inv_cumhaz_pwe(x, lambda, endpoints)
  t[curedind == 1] <- Inf
  
  # complete data observation time and event indicator
  observed.time <- pmin(t, censor.time)
  eventind      <- as.numeric(t <= censor.time)
  
  # the complete data elapsed time
  elapsed.time <- accrual.time + observed.time
  
  # if the elapsed time for a subject is larger than the administered censor time,
  # then the observed time is set to be (administered censor time - accrual time) and the event indicator is set to 0 
  idx = (elapsed.time > elapsed.time.max)
  observed.time[idx] = elapsed.time.max - (accrual.time[idx])
  eventind[idx] = 0
  
  return(list(observed.time = observed.time, eventind = eventind))
}


# Example
library(dplyr)
data.dir  <- '/proj/ibrahimlab/strapp_survival/Data'
curdata   <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
## Replace 0 times with 0.50 (half day)
curdata <- curdata %>% mutate(failtime = if_else(failtime == 0, 0.50, failtime))
range(curdata$failtime) # 0.00821 6.97604
histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))

pcure = mean(1 - histdata$fail_2yr)

n = 426
# beta, lambda and endpoints are estimated from SAS based on curdata
beta = c(-0.2063, 0.1601, -0.2112, 0.5401)
lambda = c(0.3018, 0.4234, 0.3896, 0.3765, 0.3805, 0.2002, 0.1428, 0.0695)
endpoints = c(0, 0.18754, 0.33539, 0.511975, 0.722795, 0.95551, 1.457905, 2.30801, Inf)

X <- as.matrix(cbind(curdata$trt, scale(curdata$age), curdata$sex, curdata$node_bin))
trt = X[, 1]

set.seed(1000)
res = rcurepwe(pcure = 0, X, beta, lambda, endpoints, 
               accrual.max = 52/12,
               min.follow.up = 3,
               censor.rate = log(0.95)/(-5) # solved from exp{-5x}=0.95
)

library(survival)
Y = Surv(res$observed.time, res$eventind)

par(mfrow=c(1, 2))
fit1 = survfit(Y ~ trt)
plot(fit1, lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Years", ylab = "Probability",
     main = "Relapse-Free Survival by Treatment (Simulated Data)")
# legend("topright", c("Trt", "Control"), lty = c("solid", "dashed"), col = c("black", "grey"))

fit2 = survfit(Surv(curdata$failtime, curdata$rfscens) ~ curdata$trt)
plot(fit2, lty = c("solid", "dashed"), col = c("black", "grey"), xlab = "Years", ylab = "Probability",
     main = "Relapse-Free Survival by Treatment (Real Data)")
# legend("topright", c("Trt", "Control"), lty = c("solid", "dashed"), col = c("black", "grey"))



# scaled age ~ N(0, 1)
# trt ~ Bernoulli(0.5)
# sex ~ Bernoulli(mean(histdata$sex))
# node_bin ~ Bernoulli(mean(histdata$node_bin))
# set.seed(5)
# X = cbind(rbinom(n, 1, 0.5), rnorm(n, 0, 1), rbinom(n, 1, mean(histdata$sex)), rbinom(n, 1, mean(histdata$node_bin))) # trt, scale(age), sex, node_bin
# beta = c(-0.3392, 0.1527, 0.0801, 0.6203) # beta of trt, scale(age), sex, node_bin after straPP


