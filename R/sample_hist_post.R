remove(list = ls())

library(dplyr)
library(posterior)
library(bayesplot)
library(mclust)

# ## Set cmdstanr directory (necessary for sending code to cluster)
# library(cmdstanr)
# cmdstanr::set_cmdstan_path("/nas/longleaf/home/xinxinc/.cmdstan/cmdstan-2.30.1")

## source wrappers
wrapper.dir <- '/proj/ibrahimlab/strapp_paper3/R'
source(file.path(wrapper.dir, 'wrappers.R'))

## Sampling parameters
nburnin  = 2000
nsamples = 15000
nchains  = 1

a0        <- seq(0.01, 1, by = 0.01)
num_var   <- c(1, 4) # number of covariates considered (including treatment)
grid      <- expand.grid(a0 = a0, num_var = num_var, stringsAsFactors = FALSE)
grid$seed <- sample(seq_len(1000 * nrow(grid)), nrow(grid), replace = FALSE)

## Get simulation situation based on cluster ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id) )
  id <- 1
grid.id    <- grid[id, ]
a0.id      <- grid.id$a0
num_var.id <- grid.id$num_var

# save.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_beta_4_4'
# save.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/pwe_hist_post'
save.dir <- '/proj/ibrahimlab/strapp_paper3/hist_post_samples/pwe_hist_post'

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_', 'a0_', a0.id,
                                       '_num_covariates_', num_var.id, '.rds'))

## load data
data.dir  <- '/proj/ibrahimlab/strapp_survival/Data'
histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
histdata$age_bin <- as.numeric(histdata$age >= 50)
## formula for current and historical data (including intercept)
if (num_var.id == 1) {
  fmla.hist <- fail_2yr ~ trt
} else {
  fmla.hist <- fail_2yr ~ trt + age_bin + sex + node_bin 
}

## Set seed based on task ID
set.seed(grid.id$seed)
#smpl <- hist_post(fmla.hist = fmla.hist
#                  , histdata = histdata
#                  , a0 = a0.id
#                  #, pc_upper = 0.9
#                  , isBeta = 1
#                  , pc_shape = 4
#                  , iter_warmup = nburnin, iter_sampling = nsamples
#                  , chains = 1
#                  , parallel_chains = 1)
smpl <- pwe_hist_post(fmla.hist = fmla.hist
                      , histdata = histdata
                      , a0 = a0.id
                      , iter_warmup = nburnin, iter_sampling = nsamples
                      , chains = 1
                      , parallel_chains = 1)

saveRDS(smpl, filename)

beta.strapp.vars    <- paste0('beta_strapp[', 1:num_var.id, ']')
#beta.genstrapp.vars <- paste0('beta_genstrapp[', 1:num_var.id, ']')
suppressWarnings({
  beta.strapp    <- smpl %>% as_draws_df %>% select(all_of(beta.strapp.vars)) %>% as.matrix
  #beta.genstrapp <- smpl %>% as_draws_df %>% select(all_of(beta.genstrapp.vars)) %>% as.matrix
})

## fit Gaussian mixture models 
res.strapp    <- Mclust(beta.strapp)
#res.genstrapp <- Mclust(beta.genstrapp)

res <- list(
  'id'              = id
  , 'a0'            = a0.id
  , 'num_covar'     = num_var.id
  , 'fmla.hist'     = fmla.hist
  , 'post_smpl'     = smpl
  , 'res.strapp'    = res.strapp
  #, 'res.genstrapp' = res.genstrapp
  #, 'prior_pcured'  = 'Beta(4, 4)'
)

saveRDS(res, filename)

