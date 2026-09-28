data {
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    p;             // number of covariates (EXCLUDING intercept term)
  int<lower=0,upper=1>            y0[n0];        // binary historical data
  matrix[n0,p]                    X0;            // historical data design matrix (EXCLUDING intercept term)
}

// The parameter accepted by the model
parameters {
  real                  intercept_hist;   // intercept for historical data
  vector[p]             beta_hist;        // reg. coefs for logistic reg model (excluding intercept)
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // logistic regression for historical data
  target += bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta_hist);
  // UIP
}

