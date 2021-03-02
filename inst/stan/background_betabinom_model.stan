// BETA-BINOMIAL model for NOMe-seq data from lambda DNA
// we observe that the distribution of total protected positions in free DNA (lambda DNA in our case) is overdispersed compared to binomial distribution
// here we try to model it as a beta binomial distribution 

data {
  // input data contains frequencies for occurrences for combinations of total number of positions, and number of protected positions
  int<lower=1> n_table_entries; // length of the aggregated table
  int n_total_pos[n_table_entries]; // total number of non-NA positions 
  int n_protect_pos[n_table_entries]; // total number of protected positions
  int n_frequencies[n_table_entries]; // vector containing how many combination between n_total_pos and n_protect_pos is met in the data
}


parameters {
  // parameters for beta binomial subject for inference
  real<lower=0> alpha_betabinom;
  real<lower=0> beta_betabinom;
}


model {
  
  // alpha_betabinom and beta_betabinom have uninformative non-negative priors in this model
  
  // model as a beta binom distribution
  for (n in 1:n_table_entries) {
    target += n_frequencies[n] * beta_binomial_lpmf(n_protect_pos[n] | n_total_pos[n], alpha_betabinom, beta_betabinom);
  }
}

generated quantities {
  // calculate log-likelihood function for samples
  // vector[n_table_entries] log_lik;
  // for (n in 1:n_table_entries) {
  //   log_lik[n] = n_frequencies[n] * beta_binomial_lpmf(n_protect_pos[n] | n_total_pos[n], alpha_betabinom, beta_betabinom);
  // }
  
  // mean and variance for the beta distribution
  real mean_beta;
  real var_beta;
  mean_beta = alpha_betabinom/(alpha_betabinom + beta_betabinom);
  var_beta = alpha_betabinom*beta_betabinom/(square(alpha_betabinom + beta_betabinom)*(alpha_betabinom + beta_betabinom + 1));

  
}

