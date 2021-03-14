// MIXTURE model for NOMe-seq data from lambda DNA
// we observe that the distribution of total protected positions does not match simple binomial distribution
// here we try to model it as a MIXTURE of binomial and beta binomial distribution 


// IMPORTANT NOTE! The posterior for this model is bimodal,
// because it is also possible to explain part of the data just by beta binom.
// But the second mode representing mixture of binomial and betabinomial has higher posterior
// In MAP estimate we need to use the second mode.
// Therefore, the initial point for optimisation or sampling is crucial
// source of this model is from file:
// /work2/gpeters/ozonevge/Programming/RPackages/compbio-nomeR/TEST_MODEL/01_CPP_CODE_FOR_STAN/12_mixture_model_for_lambda_NOMe.stan



data {
  // input data contains frequencies for occurrences for combinations of total number of positions, and number of protected positions
  int<lower=1> n_table_entries; // length of the aggregated table
  int n_total_pos[n_table_entries]; // total number of non-NA positions 
  int n_protect_pos[n_table_entries]; // total number of protected positions
  int n_frequencies[n_table_entries]; // vector containing how many combination between n_total_pos and n_protect_pos is met in the data
  
  // prior distribution of mixture coefficient is beta. below are params for it
  real<lower=0> alpha_distrbetaprior_mixcoeff;
  real<lower=0> beta_distrbetaprior_mixcoeff;
  
  // prior distribution of bg_protect_prob is beta. below are params for it
  real<lower=0> alpha_distrbetaprior_bg_protect;
  real<lower=0> beta_distrbetaprior_bg_protect;
  
}


parameters {
  // mixture coefficient
  real<lower=0,upper=1> theta; 
  
  // parameter for binomial component
  real<lower=0,upper=1> bg_protect_prob; 
  
  // parameters for beta binomial component
  real<lower=0> alpha_betabinom;
  real<lower=0> beta_betabinom;
}


model {
  vector[2] log_theta;
  
  // prior for mix coefficient theta
  theta ~ beta(alpha_distrbetaprior_mixcoeff,beta_distrbetaprior_mixcoeff);

  // prior for bg_protect_prob
  bg_protect_prob ~ beta(alpha_distrbetaprior_bg_protect,beta_distrbetaprior_bg_protect);

  // comment: params for beta binomial components have non-informative prior in this model version
  
  // store log of mix coeff
  log_theta[1] = log(theta);
  log_theta[2] = log(1-theta);
  // model as a mixture of binomial and beta binom distribution
  for (n in 1:n_table_entries) {
    target += n_frequencies[n] * log_sum_exp(log_theta[1] + binomial_lpmf(n_protect_pos[n] | n_total_pos[n], bg_protect_prob),
    log_theta[2] + beta_binomial_lpmf(n_protect_pos[n] | n_total_pos[n], alpha_betabinom, beta_betabinom));
  }
}

// generated quantities {
//   vector[2] log_theta;
//   // calculate log-likelihood function for samples
//   vector[n_table_entries] log_lik;
// 
//   // store log of mix coeff
//   log_theta[1] = log(theta);
//   log_theta[2] = log(1-theta);
//   // model as a mixture of binomial and beta binom distribution
//   for (n in 1:n_table_entries) {
//     log_lik[n] = n_frequencies[n] * log_sum_exp(log_theta[1] + binomial_lpmf(n_protect_pos[n] | n_total_pos[n], bg_protect_prob),
//     log_theta[2] + beta_binomial_lpmf(n_protect_pos[n] | n_total_pos[n], alpha_betabinom, beta_betabinom));
//   }
// }

