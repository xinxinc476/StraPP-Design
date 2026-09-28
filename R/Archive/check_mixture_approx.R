remove(list = ls())
library(dplyr)
library(posterior)
library(bayesplot)
library(mclust)
library(ggplot2)
library(gridExtra)

#sim.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_beta_4_4'
sim.dir <- "/pine/scr/x/i/xinxinc/strapp_paper3/curepwe_hist_post/pcured_unif_0_1"
# sim.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/pwe_hist_post'
files   <- list.files(path = sim.dir, pattern = '.rds')

sim_mixture_plot <- function(sim.res.name, method, n = 15000){
  sim.res <- readRDS(file.path(sim.dir, sim.res.name))
  smpl    <- sim.res$post_smpl
  #num_var <- length(all.vars(sim.res$fmla.hist)) - 1
  a0 = sim.res$a0
  
  if (method == "strapp"){
    beta.strapp.vars    <- paste0('beta_strapp[', 1, ']') # the coefficient for treatment
    suppressWarnings({
      beta.strapp    <- smpl %>% as_draws_df %>% select(all_of(beta.strapp.vars)) %>% as.matrix
    })
    res.strapp <- sim.res$res.strapp
    sim.strapp <- mclust::sim(res.strapp$modelName, res.strapp$parameters, n)
    df.strapp <- data.frame("sim" = sim.strapp[, 2], "post_samples" = beta.strapp[, 1])
    
    
    ncluster = length(res.strapp$parameters$pro)
    title = paste0("a0 = ", a0, ", nCluster = ", ncluster)
    
    p = ggplot(aes(x=sort(post_samples), y=sort(sim)), data=df.strapp) + 
      labs(y = "sim", x = "posterior samples", subtitle = title) +
      geom_point() +
      geom_abline(color = "red")
    
    return(p)
    
  } else{
    beta.genstrapp.vars <- paste0('beta_genstrapp[', 1, ']') # the coefficient for treatment
    suppressWarnings({
      beta.genstrapp <- smpl %>% as_draws_df %>% select(all_of(beta.genstrapp.vars)) %>% as.matrix
    })
    res.genstrapp <- sim.res$res.genstrapp
    sim.genstrapp <- mclust::sim(res.genstrapp$modelName, res.genstrapp$parameters, n)
    df.genstrapp <- data.frame("sim" = sim.genstrapp[, 2], "post_samples" = beta.genstrapp[, 1])
    
    ncluster = length(res.genstrapp$parameters$pro)
    title = paste0("a0 = ", a0, ", nCluster = ", ncluster)
    
    p = ggplot(aes(x=sort(post_samples), y=sort(sim)), data=df.genstrapp) + 
      labs(y = "sim", x = "posterior samples", subtitle = title) +
      geom_point() +
      geom_abline(color = "red")
    
    return(p) 
  }
}

set.seed(123)
id = sample(1:length(files), 9, replace = FALSE)
plots = lapply(id, function(x){sim_mixture_plot(files[x], method = "strapp")})
do.call("grid.arrange", c(plots, ncol = 3))

#plots2 = lapply(id, function(x){sim_mixture_plot(files[x], method = "genstrapp")})
#do.call("grid.arrange", c(plots2, ncol = 3))
