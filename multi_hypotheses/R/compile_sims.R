remove(list = ls())
library(dplyr)

get_power_LD <- function(sim.dir){
  files   <- list.files(path = sim.dir, pattern = '.rds')
  post_prob_betaLD = c()
  post_prob_betaHD = c()
  for (i in seq_len(length(files))) {
    simres <- readRDS(file.path(sim.dir, files[i]))$simres
    post_prob_betaLD =c(post_prob_betaLD,
                 as.numeric(unlist(simres %>% filter(variable == 'beta[1]') %>% select('post_prob')))
    )
    post_prob_betaHD =c(post_prob_betaHD,
                        as.numeric(unlist(simres %>% filter(variable == 'beta[2]') %>% select('post_prob')))
    )
  }
  
  mean(as.numeric(post_prob_betaLD > 0.975) * as.numeric(post_prob_betaHD > 0.975))
}

# n1 = 900, a0 = 0.2, strapp and pp:
get_power_LD('/pine/scr/x/i/xinxinc/strapp_paper3/multi_hypo/sims_test/strapp') # 0.5795
get_power_LD('/pine/scr/x/i/xinxinc/strapp_paper3/multi_hypo/sims_test/pp') # 0.5556
get_power_LD('/pine/scr/x/i/xinxinc/strapp_paper3/multi_hypo/sims_test/genstrapp') # 0.5662
get_power_LD('/pine/scr/x/i/xinxinc/strapp_paper3/multi_hypo/sims_test/refprior') # 0.4988


sim.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/multi_hypo/sims_test/genstrapp'
files   <- list.files(path = sim.dir, pattern = '.rds')
id      <- sapply(files, function(f){strsplit(f, '_')[[1]][2]})
id      <- as.numeric(id)
id.miss <- (400:600)[! ((400:600) %in% id) ] 

sim.dir <- '/pine/scr/x/i/xinxinc/strapp_paper3/Sims/H0/power_prior/In_progress'
files   <- list.files(path = sim.dir, pattern = '.rds')
id      <- sapply(files, function(f){strsplit(f, '_')[[1]][2]})
id      <- as.numeric(id)
for (i in 1:length(id.miss)) {
  simres <- readRDS(file.path(sim.dir, files[which(id == id.miss[i])]))$simres
  niter <- nrow(simres)/10
  print(niter)
}
