#compute average number of events for each n1 under pi_{11} and pi_{01}

remove(list = ls())

library(dplyr)

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R'
source(file.path(wrapper.dir, 'wrappers.R'))

n1.list <- seq(550, 750, by = 20)
n1.list <- c(n1.list, 1000, 1500, 2000)
n1.list <- as.integer(n1.list)
seeds_list <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/seeds_list.rds')

res.beta.strapp <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/res_beta_strapp.rds')
beta.strapp     <- res.beta.strapp$beta.strapp
X0              <- res.beta.strapp$X0
X0              <- X0[, -1] # remove the trt column
p_cured         <- res.beta.strapp$p_cured
rm(res.beta.strapp)

beta.strapp.pi11        <- c(0.85*beta.strapp[1], beta.strapp) # for pi_{11}
beta.strapp.pi01        <- c(0, beta.strapp) # for pi_{01}
beta.strapp.pi00        <- c(0, 0, beta.strapp[2:4]) # for pi_{01}
names(beta.strapp.pi00)[1:2] = names(beta.strapp.pi01)[1:2] = names(beta.strapp.pi11)[1:2] = c('trt_LD', 'trt_HD')
rm(beta.strapp)

# function to compute the average number of events for each sample size n1
# k = 1 for refprior
# k = 2 for pp and strapp
get_avg_events <- function(n1, k, beta_true){
  num_events = vector(length = length(seeds_list))
  
  for (i in seq_len(length(seeds_list))) {
    set.seed(seeds_list[i]) # set seed for generating current data
    ## Take bootstrap samples from the historical data matrix
    X.sim       <- X0[sample(1:nrow(X0), n1, replace = TRUE), ]
    curdata.sim <- rcureexp(pcure = p_cured, k = k, X = X.sim, beta = beta_true)
    num_events[i] <- sum(curdata.sim$rfscens)
  }
  mean(num_events)
}

# for refprior under pi_{11}
avg_events <- sapply(n1.list, get_avg_events, k = 1, beta_true = beta.strapp.pi11)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/refprior_pi11.rds')

# for pp/strapp under pi_{11}
avg_events <- sapply(n1.list, get_avg_events, k = 2, beta_true = beta.strapp.pi11)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/pb_pi11.rds')

# for refprior under pi_{01}
avg_events <- sapply(n1.list, get_avg_events, k = 1, beta_true = beta.strapp.pi01)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/refprior_pi01.rds')

# for pp/strapp under pi_{01}
avg_events <- sapply(n1.list, get_avg_events, k = 2, beta_true = beta.strapp.pi01)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/pb_pi01.rds')

# for refprior under pi_{00}
avg_events <- sapply(n1.list, get_avg_events, k = 1, beta_true = beta.strapp.pi00)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/refprior_pi00.rds')

# for pp/strapp under pi_{00}
avg_events <- sapply(n1.list, get_avg_events, k = 2, beta_true = beta.strapp.pi00)
names(avg_events) <- paste0('n1 = ', n1.list)
saveRDS(avg_events, file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/avg_events_per_n1/pb_pi00.rds')
