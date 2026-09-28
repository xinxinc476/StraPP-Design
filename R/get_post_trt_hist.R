# compute posterior probability of beta_trt > 0 given the historical data 
remove(list = ls())

library(dplyr)
library(posterior)
library(bayesplot)
library(cmdstanr)

## directory where stan files are located
standir <- '/proj/ibrahimlab/strapp_paper3/Stan'
hist_post_stan <- cmdstanr::cmdstan_model(file.path(standir, 'logistic_hist_post.stan'), include_paths = standir)

hist_post <- function(
  fmla.hist
  , histdata
  , ...
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
  )
  
  smpl <- hist_post_stan$sample(data = standat, ...)
  attr(smpl, 'standata')     <- standat
  attr(smpl, 'sampler.name') <- "logistic_hist_post"
  smpl
}

## load historical data
data.dir         <- '/proj/ibrahimlab/strapp_survival/Data'
histdata         <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
rm(data.dir)

fmla.hist <- fail_2yr ~ trt

set.seed(1000)
smpl <- hist_post(
  fmla.hist = fmla.hist
  , histdata = histdata
  , iter_warmup = 2000
  , iter_sampling = 25000
  , chains = 1
  , parallel_chains = 1
)

smpl$summary('beta_hist[1]', c(post_mean = mean, post_sd = sd,
                 q=~quantile(.x, probs = c(0.025, 0.975)),
                 post_prob = ~mean(.x < 0)
))
# post_prob: 0.984