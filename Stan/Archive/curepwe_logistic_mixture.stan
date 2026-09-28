
functions {
#include strapp_funs.stan
}

data {
  int<lower=0>                    n1;            // current data sample size
  int<lower=0>                    J;             // number of time intervals
  int<lower=0>                    p;             // number of covariates (no intercept term)
  vector[n1]                      y1;            // time to event current data
  matrix[n1,p]                    X1;            // current data design matrix (no intercept for current data likelihood)
  int<lower=1,upper=J>            intindx[n1];   // index giving interval into which obs i failed / was censored
  int<lower=0,upper=1>            death_ind[n1]; // event indicator (1 = event; 0 = censored)
  vector[J+1]                     breaks;        // J+1-dim interval of breaks
  real<lower=0,upper=1>           pc_upper;      // upper bound for the uniform prior of p_cured
  int<lower=0,upper=1>            isBeta;        // indicator of whether beta prior should be used for p_cured
  real<lower=0>                   pc_shape;      // shape and rate parameter for the beta prior of p_cured
  int<lower=0>                    K;             // number of components in Gaussian mixture prior
  simplex[K]                      weights;       // mixing proportions for the Gaussian mixture prior
  matrix[p, K]                    means;         // mean matrix for the Gaussian mixture prior
  array[K] cov_matrix[p]          covars;        // the kth element is the covariance matrix for the kth component of the Gaussian mixture
}

transformed data {
  real logcens[n1];
  vector[K] log_weights = log(weights);
  array[K] cov_matrix[p] precisions;
  // compute censoring indicator
  for ( i in 1:n1 )
    logcens[i] = log(1 - death_ind[i]);
  // Compute precision matrix
  for ( k in 1:K )
    precisions[k] = inverse_spd(covars[k]);
}

// The parameter accepted by the model
parameters {
  real<lower=0,upper=pc_upper> p_cured;          // proportion of cured individuals
  vector<lower=0>[J]           lambda;           // the J hazards for each interval
  vector[p]                    beta;             // reg. coefs for current data model
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[J]   log_lambda = log(lambda);
  vector[J-1] cumblhaz;
  vector[n1]  eta;
  real        lp_mix[K];
  
  // cumulative basline hazard
  cumblhaz = cumulative_sum(lambda[1:(J-1)] .* (breaks[2:J] - breaks[1:(J-1) ] ) );
    
  // half-normal prior on hazards
  target += normal_lpdf(lambda | 0, 10);
  
  // prior on p_cured
  if (isBeta == 1) {
    target += beta_lpdf(p_cured | pc_shape, pc_shape);
  } else {
    target += uniform_lpdf(p_cured | 0, pc_upper);
  }
  
  // Gaussian mixture prior: approximation to straPP/Gen-straPP prior
  for( k in 1:K ){
    // lp_mix[k] = log_weights[k] + multi_normal_lpdf(beta | means[1:p, k], covars[k]);
    lp_mix[k] = log_weights[k] + multi_normal_prec_lpdf(beta | means[, k], precisions[k]);
  }
  target += log_sum_exp(lp_mix);
  
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

