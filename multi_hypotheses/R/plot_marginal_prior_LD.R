# create density plots of marginal priors of beta_LD for straPP and pp
remove(list = ls())

library(dplyr)
library(survival)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggpubr)
library(mclust)
library(latex2exp)

# compute sd_LD
beta_HD_hat <- readRDS(file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/beta_HD_hat.rds')
sd_LD       <- 0.5*beta_HD_hat

hist.dir  <- "/proj/ibrahimlab/strapp_paper3/hist_post_samples/pwe_hist_post"
files     <- list.files(path = hist.dir, pattern = '.rds')

sample_betaLD <- function(a0, seed = 1000){
  idx        <- grepl(pattern = paste0('_a0_', a0, '_num_covariates_4'), files)
  res        <- readRDS(file.path(hist.dir, files[idx])) 
  res.strapp <- res$res.strapp
  smpl       <- res$post_smpl
  
  set.seed(seed)
  # sample betaLD for strapp
  sim.strapp    <- mclust::sim(res.strapp$modelName, res.strapp$parameters, 15000)
  betaHD.strapp <- sim.strapp[, 2]
  betaLD.strapp <- rnorm(15000, mean = 0.5*betaHD.strapp, sd = sd_LD)
  
  # sample betaLD for pp
  betaHD.pp <- smpl %>% extract_variable('beta_hist[1]')
  betaLD.pp <- rnorm(15000, mean = 0.5*betaHD.pp, sd = sd_LD)
  
  df <- data.frame(beta_LD = c(betaLD.strapp, betaLD.pp), prior = rep(c('straPP', 'PP'), each = 15000), a0 = a0)
  return(df)
}

df <- lapply(c(0.2, 0.4, 0.8), sample_betaLD)
df <- do.call(rbind, df)

df.strapp = df %>% filter(prior == 'straPP') %>% group_by(a0) %>% summarise('mean(beta_LD < 0 ) for straPP' = mean(beta_LD < 0))
df.pp = df %>% filter(prior == 'PP') %>% group_by(a0) %>% summarise('mean(beta_LD < 0 ) for PP' = mean(beta_LD < 0))
df.combined = cbind(df.strapp, df.pp[,2])
View(as.data.frame(df.combined))

appender <- function(string){TeX(paste("$a_0 = $", string))}

plt_LD = ggplot(df, aes(x = beta_LD, linetype = prior, shape = prior, color = prior))+
  labs(x = TeX("$\\beta_{LD}$"))+
  geom_density(lwd = 1.1)+
  facet_wrap(~a0, scales = "fixed", 
             labeller = as_labeller(appender, default = label_parsed))
annotate_figure(plt_LD, top = text_grob(TeX("Marginal Prior Density Plots of $\\beta_{LD}$", bold = T), 
                                        size = 14)
)
ggsave("Rplot.eps", device = "eps", path = "/proj/ibrahimlab/strapp_paper3/multi_hypotheses/Figures"
       , width = 1200, height = 480, units = "px", dpi = 300)
#width: 1200, height: 480

sample_betaHD <- function(a0, seed = 1000){
  idx        <- grepl(pattern = paste0('_a0_', a0, '_num_covariates_4'), files)
  res        <- readRDS(file.path(hist.dir, files[idx])) 
  res.strapp <- res$res.strapp
  smpl       <- res$post_smpl
  
  set.seed(seed)
  # sample betaLD for strapp
  sim.strapp    <- mclust::sim(res.strapp$modelName, res.strapp$parameters, 15000)
  betaHD.strapp <- sim.strapp[, 2]
  
  # sample betaLD for pp
  betaHD.pp <- smpl %>% extract_variable('beta_hist[1]')
  
  df <- data.frame(beta_HD = c(betaHD.strapp, betaHD.pp), prior = rep(c('straPP', 'PP'), each = 15000), a0 = a0)
  return(df)
}

df <- lapply(c(0.2, 0.4, 0.8), sample_betaHD)
df <- do.call(rbind, df)

ggplot(df, aes(x = beta_HD, linetype = prior, shape = prior, color = prior))+
  labs(x = 'beta_LD')+
  geom_density(lwd = 1.1)+
  facet_wrap(~a0, labeller = labeller(a0 = ~ paste("a0 = ", .x)) )
