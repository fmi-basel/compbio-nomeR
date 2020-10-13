#ifndef _ftp_model_loglik_start_prob_hpp_
#define _ftp_model_loglik_start_prob_hpp_


/// fully templated function for loglik function for ftp abundance
/// this version of the loglik function takes as an input vector START probabilities (not CVOER probabilties as in the previous version)
template <typename T0__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T2__, T3__>::type
  ftp_model_loglik_for_start_prob(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& ftp_start_probs,
                  const int& n_ftp,
                  const std::vector<int>& ftp_lengths,
                  const T2__& bg_protect_prob,
                  const T3__& footprint_protect_prob,
                  const std::vector<int>& spacings,
                  const int& n_spac,
                  const int& max_spacing,
                  const std::vector<std::vector<int> >& emp_joint_counts, std::ostream* pstream__)
  {
    
    // calc S_const which is for conversion ftp_cover_probs[i] = ftp_start_probs[i] * i/S_const
    T0__ S_const=0;
    for(int ftp=0;ftp < n_ftp;++ftp){
      S_const += ftp_start_probs(ftp) * (ftp + 1);
    }
    if(S_const <= 0){
      std::stringstream errmsg;
      errmsg << "ERROR! Value of S_const is negative or zero.";
      throw std::domain_error(errmsg.str());
      //Rcpp::stop("DNAbind_obj_vector::calc_theor_joint_prob: ERROR! Value of S_const is negative or zero.");
    }
    
    // calculate cover probs for ftps
    std::vector<T0__ > ftp_cover_probs;
    for(int ftp=0; ftp < n_ftp; ++ftp){
      ftp_cover_probs.push_back((ftp_start_probs(ftp) * ftp_lengths[ftp])/S_const);
    }
    
    // calc R_const which is a norm constant for conversion between cover and start probabilities (I use it in some formulas)
    T0__ R_const=0;
    for(int ftp=0;ftp < n_ftp;++ftp){
      R_const += ftp_cover_probs(ftp)/(ftp + 1);
    }
    
    if(R_const <= 0){
      std::stringstream errmsg;
      errmsg << "ERROR! Value of R_const is negative or zero.";
      throw std::domain_error(errmsg.str());
      //Rcpp::stop("DNAbind_obj_vector::calc_theor_joint_prob: ERROR! Value of R_const is negative or zero.");
    }
    
    
    
    // calculate forward partition sum and cumulative sum for sigma_cumul = \sum P(bg|d)
    std::vector<T0__ > FW(max_spacing + 1,0); // here FW[0] is initial condition for FW
    std::vector<T0__ > sigma_cumul(max_spacing + 1,0);
    FW[0] = 1;
    
    for(int d=1; d <= max_spacing; ++d){
      for(int ftp=0; ftp < n_ftp; ++ftp){
        if(d - ftp_lengths[ftp] < 0)
          break;
        FW[d] += ftp_start_probs[ftp] * FW[d - ftp_lengths[ftp]];
      }
      // calculate cumulative sum
      sigma_cumul[d] = sigma_cumul[d-1] + ftp_start_probs[0] * FW[d-1];
    }
    
    
    
    // set emission probs vectors for convenience
    T2__ beta1 = bg_protect_prob;
    T2__ beta0 = 1 - bg_protect_prob;
    T3__ alpha1 = footprint_protect_prob;
    T3__ alpha0 = 1 - footprint_protect_prob;
    
    std::vector<T3__ > alpha1_vec;
    std::vector<T3__ > alpha0_vec;
    alpha1_vec.push_back(beta1);
    alpha0_vec.push_back(beta0);
    for(int ftp=1; ftp < n_ftp; ++ftp){
      alpha1_vec.push_back(alpha1);
      alpha0_vec.push_back(alpha0);
    }
    
    
    // calculate probabilities for observing 0 and 1
    T0__ prob_tot0 = alpha0 + (beta0 - alpha0) * ftp_cover_probs[0];
    T0__ prob_tot1 = alpha1 + (beta1 - alpha1) * ftp_cover_probs[0];
    
    // add log likelihood of observing N00, N11 in the data modling using binomial distribution.
    // First row must contain data at S=0. We check that in the transformed data block of the stan file
    T0__ log_likelihood = emp_joint_counts[0][0] * log(prob_tot0) + emp_joint_counts[0][3] * log(prob_tot1);
    
    
    
    
    for(int S_idx = 1; S_idx < n_spac; ++S_idx){
      int S = spacings[S_idx];
      
      std::vector<T0__ > prob_joint(4,0); // This correspond to Spacing, P(0,0|S), P(0,1|S), P(1,0|S), P(1,1|S)
      
      // calculate sum of sigma across ftps
      T0__ sigma_sum0 = 0;
      T0__ sigma_sum1 = 0;
      
      // TODO: rewrite the for loop below to get rid of if statement
      for(int ftp=0; ftp < n_ftp; ++ftp){
        int ftp_len = ftp_lengths[ftp];
        if(S - ftp_len - 1 <= 0)
          break;
        sigma_sum0 += ftp_start_probs[ftp] * alpha0_vec[ftp] * sigma_cumul[S - ftp_len - 1];
        sigma_sum1 += ftp_start_probs[ftp] * alpha1_vec[ftp] * sigma_cumul[S - ftp_len - 1];
      }
      
      // calculate joint probabilities
      // TODO: replace R_const * ftp_start_probs[0] by ftp_cover_probs[0] after making sure that the code works as it is now
      // P(0,0|S)
      prob_joint[0] = alpha0 * alpha0 + R_const * ftp_start_probs[0] * alpha0 * (beta0 - alpha0) +
        R_const * (beta0 - alpha0) * ( (alpha0 + (beta0 - alpha0) * ftp_start_probs[0]) * sigma_cumul[S-1] - sigma_sum0);
      // P(1,1|S)
      prob_joint[3] = alpha1 * alpha1 + R_const * ftp_start_probs[0] * alpha1 * (beta1 - alpha1) +
        R_const * (beta1 - alpha1) * ( (alpha1 + (beta1 - alpha1) * ftp_start_probs[0]) * sigma_cumul[S-1] - sigma_sum1);
      // P(0,1|S)
      prob_joint[1] = (1 - prob_joint[0] - prob_joint[3])/2;
      // P(1,0|S)
      prob_joint[2] = prob_joint[1];
      
      // add to loglikelihood 
      for(int pj_ind=0; pj_ind < 4; ++pj_ind){
        log_likelihood += emp_joint_counts[S_idx][pj_ind] * log(prob_joint[pj_ind]);
      }
      
    }
    
    
    return(log_likelihood);
    
    
    
  }



#endif


