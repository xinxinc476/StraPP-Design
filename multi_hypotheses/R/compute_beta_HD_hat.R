# sd_LD = 0.5 * abs(beta_HD_hat), where beta_LD | beta_HD ~ N(0.5*beta_HD, sd_LD)
remove(list = ls())

library(dplyr)

hist.dir <- "/proj/ibrahimlab/strapp_paper3/hist_post_samples/pwe_hist_post"
files    <- list.files(path = hist.dir, pattern = '.rds')
idx      <- grepl(pattern = '_num_covariates_4', files)
files    <- files[idx]
a0_vals  <- as.numeric(sapply(files, function(x){
  strsplit(x, '_')[[1]][4]
}))

mean_mixture_list <- vector(length = length(files))
for (i in 1:length(files)) {
  res                  <- readRDS(file.path(hist.dir, files[i])) 
  res.strapp           <- res$res.strapp$parameters
  weights              <- res.strapp$pro
  means                <- res.strapp$mean
  mean_mixture_list[i] <- abs(sum(weights * means[1, ]))
}

round(quantile(mean_mixture_list), 4) # 0.1665 0.2332 0.2376 0.2395 0.2427

plot(mean_mixture_list~a0_vals, xlab = 'a0', ylab='mean of normal mixture approx', type = 'l')

mean_mixture_list[which(a0_vals == 1)] # 0.2427304
saveRDS(mean_mixture_list[which(a0_vals == 1)], 
        file = '/proj/ibrahimlab/strapp_paper3/multi_hypotheses/R/beta_HD_hat.rds')