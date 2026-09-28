
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    p;             // number of covariates (no intercept term)
  int<lower=0,upper=1>            y0[n0];        // binary historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
  real<lower=0,upper=1>           a0;            // power prior parameter
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
  vector[p] beta_strapp;
  
  // apply straPP transformation under the PWE model
  beta_strapp = invsqrt_fisher_pwe 
           * ( sqrt_fisher_logistic(X0, beta_hist, intercept_hist) * beta_hist );
}

