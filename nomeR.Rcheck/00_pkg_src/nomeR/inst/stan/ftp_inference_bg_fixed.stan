// model for footprint abundance inference in SMF data with bg_protect_prob fixed to user-defined value
/**

This model is used for bayesian inference of footprint spectrum and ftp (ftp_protect_prob) emission probability
when bg (bg_protect_prob) emission probability is fixed and defined by user.
For bg emission probability beta distribution is used as a prior.
For footprint spectrum dirichlet distribution is used as a prior.
Range of footprint lengths are selected by a user.
Log-likelihood function is calculated using a C++ function ftp_model_loglik defined in a separate hpp file in the folder inst/include

**/
functions {

	/**
   * A declaration function for calculating log-likelihood probability for a given vector of cover probs
   * @param ftp_cover_probs Vector containing cover probabilities
   * @param ftp_lengths Vector containing footprint lengths
   * @param n_ftp Length of vector ftp_cover_probs representing number of footprints,
   * @param bg_protect_prob background protect probability
   * @param ftp_protect_prob footprint protect probability
   * @param spacings vector containing spacing values (S)
   * @param n_spac length of spacings
   * @param max_spacing maximum value in vector spacings
   * @param emp_joint_counts matrix containing count data for observed combination of 00,01,10,11.
   * @return Log-likelihood value.
  */

	real ftp_model_loglik(vector ftp_cover_priors, // footprint coverages
  array[] int ftp_lengths, // 1d array which contains ftp lengths. must have equal size to ftp_cover_priors
  int n_ftp, // number of footprints, including background
  real bg_protect_prob, // emission prob for background
  real ftp_protect_prob, // emission prob for footprint
  array[] int spacings, // vector of spacings
  int n_spac, // length of the vector of spacings
  int max_spacing, // maximum spacing in spacings
  int[,] emp_joint_counts // coocurrence count table, this matrix must contain in the first row total number of 0 and 1, i.e. at S=1
  );
}


data {
  // input count coocurrnce table
  int n_spac;
  array[n_spac] int<lower=1> spacings;
  array[n_spac,4] int<lower=0> spacing_counts;

  // footprint hyperparameters
  int n_ftp;
  array[n_ftp] int ftp_lengths;
  vector[n_ftp] ftp_prior_cover; // vector of expected footprint coverages. this vector incudes 1, i.e. background and ith element represent footprint length
  real total_cnt_prior_dirich; // total count for dirichlet distribution

  // emission prob for BG for inference
  // in this version of the model this parameter is constant and defined by user
  real<lower=0,upper=1> bg_protect_prob;

  // footprint noise hyperparameters
  real ftp_protect_min;
  real ftp_protect_max;
  real ftp_protect_mean;
  real ftp_protect_totcount;
}

transformed data {
  // ftp cover dirichlet alpha hyperparams
  vector[n_ftp] ftp_prior_alphas = ftp_prior_cover * total_cnt_prior_dirich;


  // ftp noise beta distr hyperparams
  real ftp_protect_alpha = ftp_protect_totcount * ftp_protect_mean;
  real ftp_protect_beta = ftp_protect_totcount * (1 - ftp_protect_mean);


  // initialize max_spacing
  int max_spacing = max(spacings);
}



// parameters which we want to infer are footprint abundance
parameters {
  simplex[n_ftp] ftp_abundances; // footprint coverages to infer

  // emission prob for ftp for inference
  real<lower=ftp_protect_min,upper=ftp_protect_max> ftp_protect_prob;
}

model {

  // ftp abundances are modeled by dirichlet distribution
  ftp_abundances ~ dirichlet(ftp_prior_alphas);

  //prior for footprint_protect_prob
  ftp_protect_prob ~ beta(ftp_protect_alpha, ftp_protect_beta);

	// add log-likelihood
	target += ftp_model_loglik(ftp_abundances,
															ftp_lengths,
                              n_ftp,
                              bg_protect_prob,
                              ftp_protect_prob,
                              spacings,
                              n_spac,
                              max_spacing,
                              spacing_counts
  );
}




