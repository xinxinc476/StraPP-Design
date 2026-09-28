
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    n0;            // historical data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates excluding the intercept
  vector[n1]                      y1;            // time to event current data
  array[n0] int<lower=0,upper=1>  y0;            // binary historical data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  matrix[n0,p]                    X0;            // historical data design matrix (excluding intercept term)
  array[n1] int<lower=1,upper=J>  intindx;       // index giving interval into which obs i failed / was censored
  array[n1] int<lower=0,upper=1>  death_ind;     // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  real<lower=0,upper=1>           a0;            // power prior parameter
  real<lower=0,upper=1>           pc_upper;      // upper bound for the uniform prior of p_cured
  int<lower=0,upper=1>            isBeta;        // indicator of whether beta prior should be used for p_cured
  real<lower=0>                   pc_shape;      // shape and rate parameter for the beta prior of p_cured
}

transformed data {
  vector[n1] logcens;

  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
}

// The parameter accepted by the model
parameters {
  real<lower=0,upper=pc_upper> p_cured; // proportion of cured individuals
  vector<lower=0>[J]    lambda;         // the J hazards for each interval
  vector[p]             beta;           // reg. coefs
  real                  intercept_hist;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;

  // cumulative basline hazard
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );

  // half-normal(0, 10) prior on hazards
  target += normal_lpdf(lambda | 0, 10);

  // prior on p_cured
  if (isBeta == 1) {
    target += beta_lpdf(p_cured | pc_shape, pc_shape);
  } else {
    target += uniform_lpdf(p_cured | 0, pc_upper);
  }

  // power prior: logistic regression for historical data
  target += a0 * bernoulli_logit_glm_lpmf(y0 | X0, intercept_hist, beta);

  // noninformative initial prior
  target += normal_lpdf(intercept_hist | 0, 10);
  target += normal_lpdf(beta | 0, 10);

  // Current data likelihood:
  //  log_mix(p, l1, l2) = log( exp( log(p) + l1 ) + exp( log(p) + l2 ) )
  //                     = log( p * exp(l1) + (1 - p) * exp(l2) )
  //     p  = mixture prob
  //     l1 = log likelihood for first mixture  = log(cens[i]) = { log(-Inf) if event; 0 if censored }
  //     l2 = log likelihood for second mixture = log PWE likelihood
  eta = X1 * beta;
  for ( i in 1:n1 ) {
    target += log_mix(
      p_cured,
      logcens[i],
      pwe_lpdf(y1[i] | eta[i], lambda, log_lambda, breaks, intindx[i], death_ind[i], cumblhaz)
    );
  }
}


