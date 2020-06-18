    
functions{
  real ftp_model_loglik(vector ftp_cover_priors,  // here vector of priors also represent lengths, namely ith element of the vector
  // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
  int n_ftp, // number of footprints
  int[] ftp_lengths,
  real bg_protect_prob,
  real footprint_protect_prob,
  int[] spacings,// this vector must contain as first element S=1
  int n_spac,
  int max_spacing,
  int[,] emp_joint_counts //this matrix must contain in the first row total number of 0 and 1, i.e. at S=1
  );

}


data {
  int n_ftp;
  vector<lower=0>[n_ftp] ftp_prior_alphas;
  int n_spac;
  int<lower=1> spacings[n_spac];
  int<lower=0> spacing_counts[n_spac,4];
  
  // parameters for bg_protect prior
  // bg_protect_prob is const and ftp_protect variable
  real bg_protect_prob;

// parameters for ftp_protect prior. modelled using truncated beta distribution
  real ftp_protect_min;
  real ftp_protect_max;
  real ftp_protect_mean;
  real ftp_protect_var;

  
}

transformed data {
  int<lower=1> ftp_lengths[n_ftp];
  real<lower=0> ftp_protect_alpha;
  real<lower=0> ftp_protect_beta;
  int<lower=1> max_spacing;
  
  // fill ftp_lengths
  for (ftpl in 1:n_ftp){
    ftp_lengths[ftpl] = ftpl;
  }
  
  // initialize max_spacing
  max_spacing = max(spacings);
  
  // check whether first element in spacings is 1
  if(spacings[1] != 1){
    reject("Incorrect first element in spacings! First element in vector spacings (S) must be 1 which corresponds to total number of 0 and 1 in the data");
  }
  if(!(spacing_counts[1,1] > 0 && spacing_counts[1,2]==0 && spacing_counts[1,3] == 0 && spacing_counts[1,4] > 0)){
    reject("Incorrect first row in spacing_counts! It must contain N00, N01, N10, N11 at spacing S=1. Namely N00 and N11 must represent nubmer of 0 and 1 in the data and N01 and N10 must be 0.");
  }
  
  //params for FTP emis beta distr priors
  
  if(ftp_protect_var < ftp_protect_mean * (1 - ftp_protect_mean)){
    real ftp_nu = ftp_protect_mean * (1 - ftp_protect_mean)/ftp_protect_var - 1;
    ftp_protect_alpha = ftp_protect_mean * ftp_nu;
    ftp_protect_beta = (1 - ftp_protect_mean) * ftp_nu;
  } else{
    reject("Incorrect values for ftp_protect_mean and ftp_protect_var! Inequality ftp_protect_var < ftp_protect_mean * (1 - ftp_protect_mean) must hold.");
  }
  
}

// parameters which we want to infer are footprint abundance
parameters {
  simplex[n_ftp] ftp_abundances;
  
  // emission prob for ftp
  real<lower=ftp_protect_min,upper=ftp_protect_max> footprint_protect_prob;
}

model {
  // define prior for abundances
  ftp_abundances ~ dirichlet(ftp_prior_alphas);

  //define prior for footprint_protect_prob
  footprint_protect_prob ~ beta(ftp_protect_alpha, ftp_protect_beta);
  
  // add log likelihood
  target += ftp_model_loglik(ftp_abundances,
                              n_ftp,
                              ftp_lengths,
                              bg_protect_prob,
                              footprint_protect_prob,
                              spacings,
                              n_spac,
                              max_spacing,
                              spacing_counts
  );
}
