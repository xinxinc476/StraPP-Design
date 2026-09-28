remove(list = ls())

library(dplyr)
library(survival)
library(posterior)
library(bayesplot)

# ## Set cmdstanr directory (necessary for sending code to cluster)
# library(cmdstanr)
# cmdstanr::set_cmdstan_path("/nas/longleaf/home/xinxinc/.cmdstan/cmdstan-2.30.1")

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_paper3/R'
source(file.path(wrapper.dir, 'wrappers.R'))

## load data
data.dir  <- '/proj/ibrahimlab/strapp_survival/Data'
histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
curdata   <- readRDS(file.path(data.dir, 'data_cur_1690.rds'))
## Replace 0 times with 0.50 (half day)
curdata <- curdata %>% mutate(failtime = if_else(failtime == 0, 0.50, failtime))
nbreaks = 5
probs   = 1:nbreaks / nbreaks
breaks  = curdata %>%
  filter(rfscens == 1) %>%
  summarize(quant = quantile(failtime, probs = probs)) %>%
  unlist
breaks = c(0, breaks)
breaks[length(breaks)] = max(10000, 1000 * breaks[length(breaks)])

# sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_beta_4_4"
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/pwe_hist_post"
files   <- list.files(path = sim.dir, pattern = '.rds')

# ncomponents = vector(length = length(files))
# for(i in 1:length(files)){
#  res <- readRDS(file.path(sim.dir, files[i]))
#  res.strapp <- res$res.strapp$parameters
#  weights <- res.strapp$pro
#  ncomponents[i] = length(weights)
#}
# table(ncomponents)
# all > 1, 49 have 2 components, 85 have 3 components, 30 have 4 components, 20 have 5 components 

#res <- readRDS(file.path(sim.dir, files[4])) # have 4 covariates
res <- readRDS(file.path(sim.dir, files[2])) # only have trt as the covariate
fmla.hist <- res$fmla.hist
if (res$num_covar == 1) {
  fmla.cur  <- Surv(failtime, rfscens) ~ trt
} else {
  fmla.cur  <- Surv(failtime, rfscens) ~ trt + scale(age) + sex + node_bin
}

res.strapp <- res$res.strapp$parameters
weights <- res.strapp$pro
means <- res.strapp$mean
covars <- (res.strapp$variance)[[4]]

nburnin  = 2000
nsamples = 15000
smpl <- curepwe_mix(fmla.cur = fmla.cur
                  , data = curdata
                  , breaks = breaks
                  , pc_upper = 1
                  , isBeta = 1
                  , pc_shape = 4
                  , weights = weights
                  , means = means
                  , covars = covars
                  , iter_warmup = nburnin, iter_sampling = nsamples
                  , chains = 1
                  , parallel_chains = 1)

smpl.strapp <- curepwe_mix_strapp(fmla.cur = fmla.cur
                    , data = curdata
                    , breaks = breaks
                    , pc_upper = 1
                    , isBeta = 1
                    , pc_shape = 4
                    , weights = weights
                    , means = means
                    , covars = covars
                    , iter_warmup = nburnin, iter_sampling = nsamples
                    , chains = 1
                    , parallel_chains = 1)

smpl.genstrapp <- curepwe_mix_genstrapp(fmla.hist = fmla.hist
                                  , fmla.cur = fmla.cur
                                  , histdata = histdata
                                  , data = curdata
                                  , breaks = breaks
                                  , pc_upper = 1
                                  , isBeta = 1
                                  , pc_shape = 4
                                  , weights = weights
                                  , means = means
                                  , covars = covars
                                  , iter_warmup = nburnin, iter_sampling = nsamples
                                  , chains = 1
                                  , parallel_chains = 1)
