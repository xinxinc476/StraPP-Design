## simulation scenario: under H0, gen-strapp
remove(list = ls())

library(dplyr)
library(survival)
library(posterior)
library(bayesplot)

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_paper3/R'
source(file.path(wrapper.dir, 'wrappers.R'))

## Sampling parameters
nburnin  = 2000
nsamples = 25000
nchains  = 1

n1 <- seq(550, 750, by = 20)
a0 <- seq(0.2, 0.8, by = 0.2)

grid <- expand.grid(
  n1 = n1,
  a0 = a0,
  stringsAsFactors = FALSE
)

ndatasets <- 10000  ## total number of data sets
each.cl   <- 100 ## how many data sets to run on a single node
ncl       <- ceiling(ndatasets / each.cl) ## how many nodes per data set

## Repeat each row of the grid ncl times (want nrow(grid) <= 800 ideally)
grid <- grid[rep(1:nrow(grid), each = ncl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1)

# generate seeds for each dataset
seeds_list <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/seeds_list.rds')

## Get simulation situation based on cluster ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id) )
  id <- 1
grid.id  <- grid[id, ]
n1.id    <- grid.id$n1
a0.id    <- grid.id$a0
seeds.id <- seeds_list[grid.id$start:grid.id$end]

## Obtain file name 
save.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/Sims/H0/genstrapp'
filename <- file.path(save.dir, paste0('id_', id, '_', 'n1_', n1.id, '_a0_', a0.id, '_', grid.id$start, '-', grid.id$end, '.rds'))

## load mean and covariance components for the Gaussian mixture approximation
hist.dir  <- "/proj/ibrahimlab/strapp_paper3/hist_post_samples/pwe_hist_post"
files     <- list.files(path = hist.dir, pattern = '.rds')
idx       <- grepl(pattern = paste0('_a0_', a0.id, '_num_covariates_4'), files)
res       <- readRDS(file.path(hist.dir, files[idx])) 
fmla.hist <- res$fmla.hist

res.strapp <- res$res.strapp$parameters
weights <- res.strapp$pro
means <- res.strapp$mean
covars <- (res.strapp$variance)[[4]]

# load beta.strapp and X0 (data matrix from historical data)
res.beta.strapp <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/res_beta_strapp.rds')
beta.strapp     <- res.beta.strapp$beta.strapp
beta.strapp[1]  <- 0 # for null hypothesis
X0              <- res.beta.strapp$X0
p_cured         <- res.beta.strapp$p_cured
rm(res.beta.strapp)

## formula for current data 
fmla.cur  <- Surv(failtime, rfscens) ~ trt + age_bin + sex + node_bin

## load historical data
data.dir         <- '/proj/ibrahimlab/strapp_survival/Data'
histdata         <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
histdata$age_bin <- as.numeric(histdata$age >= 50)
rm(data.dir)

## function to obtain post_mean, post_sd, CI, P(trteff < 0) = mean(trteff < 0)
getStats <- function(fit) {
  var_names <- names(fit$draws(format = 'draws_df'))
  var_names <- var_names[var_names %in% c('p_cured', paste0('beta[', 1:4, ']'), paste0('lambda[', 1:5, ']'))]
  res = lapply(var_names, function(x){
    fit$summary(x, c(post_mean = mean, post_sd = sd,
                     q=~quantile(.x, probs = c(0.025, 0.975)),
                     post_prob = ~mean(.x < 0)
    ))
  })
  res       <- do.call(rbind, res)
  return(res)
}

for(i in seq_len(each.cl)){
  set.seed(seeds.id[i]) # set seed for generating current data
  ## Take bootstrap samples from the historical data matrix
  X.sim       <- X0[sample(1:nrow(X0), n1.id, replace = TRUE), ]
  curdata.sim <- rcureexp(pcure = p_cured, X = X.sim, beta = beta.strapp)
  
  nbreaks = 5
  probs   = 1:nbreaks / nbreaks
  breaks  = curdata.sim %>%
    filter(rfscens == 1) %>%
    summarize(quant = quantile(failtime, probs = probs)) %>%
    unlist
  breaks = c(0, breaks)
  breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])
  
  smpl <- curepwe_mix_genstrapp(fmla.hist = fmla.hist
                                , fmla.cur = fmla.cur
                                , histdata = histdata
                                , data = curdata.sim
                                , breaks = breaks
                                , pc_upper = 1
                                , isBeta = 1
                                , pc_shape = 4
                                , weights = weights
                                , means = means
                                , covars = covars
                                , iter_warmup = nburnin, iter_sampling = nsamples
                                , chains = 1
                                , parallel_chains = 1
                                , adapt_delta = 0.99
                                , max_treedepth = 12
  ) 
  
  simres.i <- getStats(smpl)
  
  if(i == 1) {
    simres <- simres.i
  } else{
    simres <- rbind(simres, simres.i)
  }
  
  if(i %% 20 == 0) {
    saveRDS(list(
      'id'       = id
      , 'n1'     = n1.id
      , 'a0'     = a0.id
      , 'simres' = simres
      , 'prior'  = 'genstrapp'
      #, 'hypothesis' = 'alternative'
      , 'hypothesis' = 'null'
    ), file = file.path(save.dir,
                        paste0('In_progress/', 'id_', id, '_', 'n1_', n1.id, '_a0_', a0.id, '_', grid.id$start, '-', grid.id$end, '.rds')))
  }
  print(paste0("#################### Completed iteration ", i, " ######################"))
}

res <- list(
  'id'       = id
  , 'n1'     = n1.id
  , 'a0'     = a0.id
  , 'simres' = simres
  , 'prior'  = 'genstrapp'
  #, 'hypothesis' = 'alternative'
  , 'hypothesis' = 'null'
)

saveRDS(res, filename)



