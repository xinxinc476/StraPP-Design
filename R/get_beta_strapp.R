remove(list = ls())

library(dplyr)

# function to compute the inverse square root of a symmetric matrix
inv_sqrt_symm <- function(mtx){
  res       <- eigen(mtx)
  eigen_val <- res$values
  V         <- res$vectors
  return( V %*% diag(1 / sqrt(eigen_val), nrow = dim(mtx)[1]) %*% t(V) )
}

# function to compute the square root of a symmetric matrix
sqrt_symm <- function(mtx){
  res       <- eigen(mtx)
  eigen_val <- res$values
  V         <- res$vectors
  return( V %*% diag(sqrt(eigen_val), nrow = dim(mtx)[1]) %*% t(V) )
}

# function to return square root fisher information matrix dropping intercepts
# beta are the MLEs from logistic regression (Not including the intercept)
sqrt_fisher_logistic <- function(X, beta, intercept) {
  W = as.numeric(binomial(link = 'logit')$linkinv(intercept + X %*% beta))
  W = W * (1 - W)
  W_sqrt = diag(sqrt(W), nrow = length(W))
  return(sqrt_symm(crossprod(W_sqrt %*% X)))
}

data.dir         <- '/proj/ibrahimlab/strapp_survival/Data'
histdata         <- readRDS(file.path(data.dir, 'data_hist_1684.rds'))
histdata$age_bin <- as.numeric(histdata$age >= 50)
p_cured          <- mean(1 - histdata$fail_2yr)

fmla.hist <- fail_2yr ~ trt + age_bin + sex + node_bin 
# fmla.hist <- fail_2yr ~ trt
fit.glm   <- glm(fmla.hist, family = binomial(link = 'logit'), data = histdata)
beta.glm  <- coef(fit.glm) # MLEs from logistic regression
intercept <- beta.glm[1]
beta.glm  <- beta.glm[-1]

X0        <- model.matrix(fmla.hist, histdata)
## Make sure no design matrices have intercepts
if ( '(Intercept)' %in% colnames(X0) ) {
  X0 <- X0[, -1, drop = FALSE]
}
invsqrt_fisher_pwe = inv_sqrt_symm(crossprod(X0))

# apply straPP transformation
# under the alternative
beta.strapp <- (invsqrt_fisher_pwe %*% sqrt_fisher_logistic(X0, beta.glm, intercept) %*% beta.glm) / sqrt(1 - p_cured)
beta.strapp <- as.numeric(beta.strapp)
names(beta.strapp) <- names(beta.glm)
round(beta.strapp, 4) # -0.3217   0.3711   0.0794   0.6436

# save p_cured, beta.strapp and X0
saveRDS(list(p_cured = p_cured, beta.strapp = beta.strapp, X0 = X0), file = '/proj/ibrahimlab/strapp_paper3/R/res_beta_strapp.rds')

beta.strapp.null <- beta.strapp
beta.strapp.null[1] <- 0
round(beta.strapp.null, 4)
