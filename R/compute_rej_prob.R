# plot histogram of posterior probabilities for n1 = 10000 under H0, refprior
remove(list = ls())
library(dplyr)

#res.beta.strapp <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/R/res_beta_strapp.rds')
#beta.strapp     <- res.beta.strapp$beta.strapp
#var_names = c('beta[1]', 'beta[2]', 'beta[3]', 'beta[4]') # covariate names in simres

get_rej_prob <- function(sim.dir){
  files   <- list.files(path = sim.dir, pattern = '.rds')
  post_prob = c()
  for (i in seq_len(length(files))) {
    simres <- readRDS(file.path(sim.dir, files[i]))$simres
    post_prob =c(post_prob,
                 as.numeric(unlist(simres %>% filter(variable == 'beta[1]') %>% select('post_prob')))
    )
  }
  
  return(post_prob)
}

sim.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/sims_test/refprior/check_type_I_err/n1_10000_nMCMC_25000_trt_only'
post_prob <- get_rej_prob(sim.dir = sim.dir)
mean(post_prob > (1 - 0.025)) 

library(ggplot2)
df <- data.frame('post_prob' = post_prob)

ggplot(df, aes(x=post_prob)) + 
  labs(x = 'posterior probability of beta_trt < 0')+
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, color = 'red')
