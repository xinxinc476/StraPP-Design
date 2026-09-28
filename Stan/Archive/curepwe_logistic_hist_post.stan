
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    p;             // number of covariates (no intercept term)
  int<lower=0,upper=1>            y0[n0];        // binary historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  real<lower=0,upper=1>           a0;            // power prior parameter
  real<lower=0,upper=1>           pc_upper;      // upper bound for the uniform prior of p_cured
  int<lower=0,upper=1>            isBeta;        // indicator of whether beta prior should be used for p_cured
  real<lower=0>                   pc_shape;      // shape and rate parameter for the beta prior of p_cured
}

transformed data {
  matrix[p,p] invsqrt_fisher_pwe = matrix_sqrt( inverse_spd( crossprod(X0) ) );
}

// The parameter accepted by the model
parameters {
  real                  intercept_hist;  // intercept for historical data
  vector[p]             beta_hist;       // reg. coefs for logistic reg model (includes intercept)
}

model {
  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta_hist | 0, 10);
  
  // power prior: logistic regression for historical data
  if (a0 > 0) {
    target += a0 * bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta_hist);
  }
}

generated quantities {
  real p_cured;
  vector[p] beta_strapp;
  vector[p] beta_genstrapp;
  vector[p] c0;
  real<lower=0> c0_sd = abs(std_normal_rng());
  
  if (isBeta == 1) {
    p_cured = beta_rng(pc_shape, pc_shape);
  } else {
    p_cured = uniform_rng(0, pc_upper);
  }
  
  // apply straPP transformation
  beta_strapp =   
        inv_sqrt(1 - p_cured) * invsqrt_fisher_pwe 
      * ( sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist );
      
  // apply Gen-straPP transformation
  for (i in 1:p) {
    c0[i] = normal_rng(0, c0_sd); // gen strapp location parameter
  }
  beta_genstrapp =   
        inv_sqrt(1 - p_cured) * invsqrt_fisher_pwe 
      * ( sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist - c0 );
}


