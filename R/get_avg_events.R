#compute average number of events for each n1

remove(list = ls())

library(dplyr)

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_paper3/R'
source(file.path(wrapper.dir, 'wrappers.R'))

seeds_list <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/seeds_list.rds')

res.beta.strapp <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/res_beta_strapp.rds')
beta.strapp     <- res.beta.strapp$beta.strapp
X0              <- res.beta.strapp$X0
p_cured         <- res.beta.strapp$p_cured
rm(res.beta.strapp)

# function to compute the average number of events for each sample size n1
get_avg_events <- function(n1, beta.true){
  num_events = vector(length = length(seeds_list))
  
  for (i in seq_len(length(seeds_list))) {
    set.seed(seeds_list[i]) # set seed for generating current data
    ## Take bootstrap samples from the historical data matrix
    X.sim       <- X0[sample(1:nrow(X0), n1, replace = TRUE), ]
    curdata.sim <- rcureexp(pcure = p_cured, X = X.sim, beta = beta.true)
    num_events[i] <- sum(curdata.sim$rfscens)
  }
  mean(num_events)
}

n1.list    <- seq(550, 750, by = 20)
avg_events <- sapply(n1.list, get_avg_events, beta.true = beta.strapp)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/R/avg_events_per_n1_H1.rds')

beta.strapp[1] <- 0
avg_events <- sapply(n1.list, get_avg_events, beta.true = beta.strapp)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/R/avg_events_per_n1_H0.rds')
