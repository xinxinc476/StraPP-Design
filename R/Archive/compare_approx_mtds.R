remove(list = ls())

library(dplyr)
library(posterior)
library(bayesplot)
library(mclust)
library(ggplot2)
library(gridExtra)

################################## Compare straPP approximation ###########################################

## directory of results from sampling from curepwe_logistic_hist_post.stan
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_beta_4_4"
files   <- list.files(path = sim.dir, pattern = '.rds')
res1    <- readRDS(file.path(sim.dir, files[2]))
smpl1   <- res1$post_smpl

## directory of results from sampling from pwe_logistic_hist_post.stan
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/pwe_hist_post"
files   <- list.files(path = sim.dir, pattern = '.rds')
res2    <- readRDS(file.path(sim.dir, files[2]))
smpl2   <- res2$post_smpl

beta.strapp.vars    <- paste0('beta_strapp[', 1, ']') # the coefficient for treatment
suppressWarnings({
  beta.strapp    <- smpl1 %>% as_draws_df %>% select(all_of(beta.strapp.vars)) %>% as.matrix
  beta.strapp2   <- smpl2 %>% as_draws_df %>% select(all_of(beta.strapp.vars)) %>% as.matrix
})

res.strapp1 <- res1$res.strapp
res.strapp2 <- res2$res.strapp
params      <- res.strapp2$parameters
means       <- params$mean
covars      <- params$variance$sigmasq

# draw samples using the previous approximation method
set.seed(1)
sim.strapp1 <- mclust::sim(res.strapp1$modelName, res.strapp1$parameters, 15000)
sim.strapp1 <- sim.strapp1[, 2]

# draw samples using the new approximation method
set.seed(1)
p_cured     <- rbeta(15000, 4, 4)
sim.strapp2 <- sapply(1:15000, function(i){
  params$mean             <- means / sqrt(1 - p_cured[i])
  params$variance$sigmasq <- covars / (1 - p_cured[i])
  res                     <- mclust::sim(res.strapp2$modelName, params, 1)
  as.numeric(res[, 2])
})
post_samples2 <- (beta.strapp2[,1]) / sqrt(1 - p_cured)

df.strapp <- data.frame("sim1" = sim.strapp1, "sim2" = sim.strapp2, 
                        "post_samples" = beta.strapp[, 1], "post_samples2" = post_samples2)

p1 = ggplot(aes(y=sort(post_samples), x=sort(sim1)), data=df.strapp) + 
  labs(x = "Previous Approximation", y = "MCMC samples",
       subtitle = paste0("a0 = ", res1$a0, ", number of covariates = ", res1$num_covar)) +
  geom_point() +
  geom_abline(color = "red")

p2 = ggplot(aes(y=sort(post_samples), x=sort(sim2)), data=df.strapp) + 
  labs(x = "New Approximation", y = "MCMC samples",
       subtitle = paste0("a0 = ", res2$a0, ", number of covariates = ", res2$num_covar)) +
  geom_point() +
  geom_abline(color = "red")

grid.arrange(p1, p2, ncol = 2)
###########################################################################################################


################################## Compare Gen-straPP approximation ###########################################
remove(list = ls())

## directory of results from sampling from curepwe_logistic_hist_post.stan
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_beta_4_4"
files   <- list.files(path = sim.dir, pattern = '.rds')
res1    <- readRDS(file.path(sim.dir, files[2]))
smpl1   <- res1$post_smpl

## directory of results from sampling from pwe_logistic_hist_post.stan
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/pwe_hist_post"
files   <- list.files(path = sim.dir, pattern = '.rds')
res2    <- readRDS(file.path(sim.dir, files[2]))
smpl2   <- res2$post_smpl

beta.genstrapp.vars    <- paste0('beta_genstrapp[', 1, ']') # the coefficient for treatment
suppressWarnings({
  beta.genstrapp    <- smpl1 %>% as_draws_df %>% select(all_of(beta.genstrapp.vars)) %>% as.matrix
})

res.genstrapp1 <- res1$res.genstrapp
res.strapp2 <- res2$res.strapp
params      <- res.strapp2$parameters
means       <- params$mean
covars      <- params$variance$sigmasq

# draw samples using the previous approximation method
set.seed(1000)
sim.genstrapp1 <- mclust::sim(res.genstrapp1$modelName, res.genstrapp1$parameters, 15000)
sim.genstrapp1 <- sim.genstrapp1[, 2]


# load historical data
data.dir  <- '/proj/ibrahimlab/strapp_survival/Data'
histdata  <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
fmla.hist <- res2$fmla.hist
X0        <- model.matrix(fmla.hist, histdata)
## Make sure no design matrices have intercepts
if ( '(Intercept)' %in% colnames(X0) ) {
  X0 <- X0[, -1, drop = FALSE]
}
invsqrt_fisher_pwe = as.numeric(1 / sqrt( t(X0) %*% X0 ) )

# draw samples using the new approximation method
set.seed(1000)
p_cured     <- rbeta(15000, 4, 4)
sim.genstrapp2 <- sapply(1:15000, function(i){
  c0_sd <- rnorm(1)
  while (c0_sd <= 0) {
    c0_sd <- rnorm(1)
  }
  c0 <- rnorm(1, mean = 0, sd = c0_sd)
  params$mean             <- (means - invsqrt_fisher_pwe * c0) / sqrt(1 - p_cured[i])
  params$variance$sigmasq <- covars / (1 - p_cured[i])
  res                     <- mclust::sim(res.strapp2$modelName, params, 1)
  as.numeric(res[, 2])
})

df.genstrapp <- data.frame("sim1" = sim.genstrapp1, "sim2" = sim.genstrapp2, 
                        "post_samples" = beta.genstrapp[, 1])

p1 = ggplot(aes(y=sort(post_samples), x=sort(sim1)), data=df.genstrapp) + 
  labs(x = "Previous Approximation", y = "MCMC samples",
       subtitle = paste0("a0 = ", res1$a0, ", number of covariates = ", res1$num_covar)) +
  geom_point() +
  geom_abline(color = "red")

p2 = ggplot(aes(y=sort(post_samples), x=sort(sim2)), data=df.genstrapp) + 
  labs(x = "New Approximation", y = "MCMC samples",
       subtitle = paste0("a0 = ", res2$a0, ", number of covariates = ", res2$num_covar)) +
  geom_point() +
  geom_abline(color = "red")

grid.arrange(p1, p2, ncol = 2)


