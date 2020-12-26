/**
This model for footprint inference solves non-identifiablility problem 
for bg_protect_prob, ftp_protect_prob and ftp_cover_prob[1] by fixing bg_protect_prob to predefined value.


Original source of this file:
/work2/gpeters/ozonevge/Programming/RPackages/compbio-nomeR/TEST_MODEL/01_CPP_CODE_FOR_STAN/13_ftp_inference_background_fixed.stan

**/


functions{
  /**
   * function for converting vector of ftp start probabilities to cover probabilities
   * @param ftp_start_probs Vector containing start probabilities
   * @param n_ftp Length of vector ftp_start_probs representing number of footprints,
   * @param ftp_lengths Vector containing footprint lengths
   * @return vector of cover probabilities.
  */
  vector start_probs_2_cover_probs(vector ftp_start_probs,
                                   int n_ftp,
                                   int[] ftp_lengths){
    vector[n_ftp] ftp_cover_probs;
    
    for(ftp in 1:n_ftp){
      ftp_cover_probs[ftp] = ftp_start_probs[ftp] * ftp_lengths[ftp];
    }
    ftp_cover_probs = ftp_cover_probs/sum(ftp_cover_probs);
    return(ftp_cover_probs);
  }
  
  
  /**
   * function for converting vector of ftp cover probabilities to start probabilities
   * @param ftp_cover_probs Vector containing cover probabilities
   * @param n_ftp Length of vector ftp_cover_probs representing number of footprints,
   * @param ftp_lengths Vector containing footprint lengths
   * @return Vector of start probabilities.
  */
  vector cover_probs_2_start_probs(vector ftp_cover_probs,
                                   int n_ftp,
                                   int[] ftp_lengths){
    vector[n_ftp] ftp_start_probs;
    
    for(ftp in 1:n_ftp){
      ftp_start_probs[ftp] = ftp_cover_probs[ftp]/ftp_lengths[ftp];
    }
    ftp_start_probs = ftp_start_probs/sum(ftp_start_probs);
    return(ftp_start_probs);
  }
  
  /**
   * A declaration function for calculating log-likelihood probability for a given vector of cover probs
   * @param ftp_cover_probs Vector containing cover probabilities
   * @param n_ftp Length of vector ftp_cover_probs representing number of footprints,
   * @param ftp_lengths Vector containing footprint lengths
   * @param bg_protect_prob background protect probability
   * @param footprint_protect_prob footprint protect probability
   * @param spacings vector containing spacing values (S)
   * @param n_spac length of spacings
   * @param max_spacing maximum value in vector spacings
   * @param emp_joint_counts matrix containing count data for observed combination of 00,01,10,11.
   * @return Log-likelihood value.
  */
  
  real ftp_model_loglik(vector ftp_cover_probs,  // here vector of priors also represent lengths, namely ith element of the vector
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
  
  
  // emission prob for BG for inference
  // in this version of the model this parameter is constant
  real<lower=0,upper=1> bg_protect_prob;
  
  // parameters for ftp_protect prior. modelled using truncated beta distribution
  real ftp_protect_min;
  real ftp_protect_max;
  real<lower=0> ftp_protect_alpha;
  real<lower=0> ftp_protect_beta;

  
}

transformed data {
  int<lower=1> ftp_lengths[n_ftp];
  
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
}

// parameters which we want to infer are footprint abundance
parameters {
  // footprint cover probs for inference
  simplex[n_ftp] ftp_cover_probs;
  
  // emission prob for ftp for inference
  real<lower=ftp_protect_min,upper=ftp_protect_max> footprint_protect_prob;
}

model {
  real total_prob_of_protected_pos;
  // prior for abundances
  ftp_cover_probs ~ dirichlet(ftp_prior_alphas);

  //prior for footprint_protect_prob
  footprint_protect_prob ~ beta(ftp_protect_alpha, ftp_protect_beta);
  
  // add to loglik the term which realtes occurences of 1, bg_cover_prob, bg_protect_prob, ftp_protect_prob
  total_prob_of_protected_pos = bg_protect_prob * ftp_cover_probs[1] + footprint_protect_prob * (1 - ftp_cover_probs[1]);
  target += binomial_lpmf(spacing_counts[1,4] | spacing_counts[1,1] + spacing_counts[1,4], total_prob_of_protected_pos);
  
  
  // add log likelihood the window crosscorel
  target += ftp_model_loglik(ftp_cover_probs,
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

generated quantities {
  vector[n_ftp] ftp_start_probs = cover_probs_2_start_probs(ftp_cover_probs,
                                                            n_ftp,
                                                            ftp_lengths);
}
