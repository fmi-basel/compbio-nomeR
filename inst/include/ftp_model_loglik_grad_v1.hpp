#ifndef _ftp_model_loglik_and_grad_hpp_
#define _ftp_model_loglik_and_grad_hpp_


#include <stan/math/rev/core.hpp>
#include <stan/math/rev.hpp>
#include <omp.h>



// function calculating lolikelihood function as well as gradient w.r.t to footprint abundances.
// NOTE: no gradient for ftp_emiss_prob



namespace ftpmodfunctions {

double kronecker_delta(int i, int j){
  return(i == j ? 1 : 0);
}

double ftp_model_loglik_grad(const Eigen::Matrix<double, -1, 1>& ftp_cover_priors,
                             const int& n_ftp, // number of footprint. must equal to ftp_cover_priors.size()
                             const double& bg_protect_prob, // emission prob for bg
                             const double& footprint_protect_prob, // emission pros for footprints
                             const std::vector<int >& spacings, // vector containing of spacing values
                             const int& n_spac, // lengths of spacings vector. must be equal to spacings.size()
                             const int& max_spacing, // maximum spacing in vector spacings. must be equal to max(spacings)
                             const std::vector<std::vector<int > >& emp_joint_counts, // arrary (vecotr of vector) containing data for observing 00,01,10,11 at spacings
                             
                             
                             std::vector<double >& log_likelihood_gradient, // mutable vector for storing gradient, must be a vector of size n_ftp
                             //Eigen::VectorXd& log_likelihood_gradient,
                             // columns must be N00, N01, N10, N11
                             std::ostream* pstream__
){
  
  
  // calc R_const which is a norm constant for conversion between cover and start probabilities
  std::vector<int > ftp_lengths;
  double R_const=0;
  for(int ftp=0;ftp < n_ftp;++ftp){
    ftp_lengths.push_back(ftp + 1);
    R_const += ftp_cover_priors(ftp)/(ftp + 1);
  }
  
  if(R_const <= 0){
    std::stringstream errmsg;
    errmsg << "ERROR! Value of R_const is negative or zero.";
    throw std::domain_error(errmsg.str());
    //Rcpp::stop("DNAbind_obj_vector::calc_theor_joint_prob: ERROR! Value of R_const is negative or zero.");
  }
  
  
  // calculate start priors for ftps
  std::vector<double > ftp_start_priors;
  for(int ftp=0; ftp < n_ftp;++ftp){
    ftp_start_priors.push_back((ftp_cover_priors(ftp)/ftp_lengths[ftp])/R_const);
  }
  //std::cout<<"Start priors are calculated"<<std::endl;
  // calculate forward partition sum and cumulative sum for sigma_cumul = \sum P(bg|d)
  std::vector<double > FW(max_spacing + 1,0); // here FW[0] is initial condition for FW
  std::vector<double > sigma_cumul(max_spacing + 1,0);
  // calculate partition sum and sigma_cumul
  FW[0] = 1;
  //ProbBg[0] = 0;
  //sigma_cumul[0] = 0;
  //std::cout<<"Start calc FW"<<std::endl;
  for(int d=1; d <= max_spacing; ++d){
    //FW[d] = 0;
    for(int ftp=0; ftp < n_ftp; ++ftp){
      if(d - ftp_lengths[ftp] < 0)
        break;
      
      FW[d] += ftp_start_priors[ftp] * FW[d - ftp_lengths[ftp]];
    }
    // calculate cumulative sum
    sigma_cumul[d] = sigma_cumul[d-1] + ftp_start_priors[0] * FW[d-1];
  }
  
  // set emission probs vectors for convenience
  double beta1 = bg_protect_prob;
  double beta0 = 1 - bg_protect_prob;
  double alpha1 = footprint_protect_prob;
  double alpha0 = 1 - footprint_protect_prob;
  std::vector<double > alpha1_vec;
  std::vector<double > alpha0_vec;
  alpha1_vec.push_back(beta1);
  alpha0_vec.push_back(beta0);
  for(int ftp=1; ftp < n_ftp; ++ftp){
    alpha1_vec.push_back(alpha1);
    alpha0_vec.push_back(alpha0);
  }
  
  // calculate joint theoretical probabilities
  std::vector<std::vector<double > > Pjoint_array;
  //Pjoint_array.clear()
  
  
  /// calculate loglikelihood
  double log_likelihood = 0;
  
  for(int S_idx = 0; S_idx < n_spac; ++S_idx){
    int S = spacings[S_idx];
    
    std::vector<double > prob_joint(4,0); // This correspond to Spacing, P(0,0|S), P(0,1|S), P(1,0|S), P(1,1|S)
    double sigma_sum0 = 0;
    double sigma_sum1 = 0;
    
    
    for(int ftp=0; ftp < n_ftp; ++ftp){
      int ftp_len = ftp_lengths[ftp];
      if(S - ftp_len - 1 <= 0)
        break;
      sigma_sum0 += ftp_start_priors[ftp] * alpha0_vec[ftp] * sigma_cumul[S - ftp_len - 1];
      sigma_sum1 += ftp_start_priors[ftp] * alpha1_vec[ftp] * sigma_cumul[S - ftp_len - 1];
      
    }
    
    //prob_joint[0] = S;
    // calc probs
    // P(0,0|S)
    prob_joint[0] = alpha0 * alpha0 + R_const * ftp_start_priors[0] * alpha0 * (beta0 - alpha0) +
      R_const * (beta0 - alpha0) * ( (alpha0 + (beta0 - alpha0) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum0);
    
    // P(1,1|S)
    prob_joint[3] = alpha1 * alpha1 + R_const * ftp_start_priors[0] * alpha1 * (beta1 - alpha1) +
      R_const * (beta1 - alpha1) * ( (alpha1 + (beta1 - alpha1) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum1);
    
    // P(0,1|S)
    
    prob_joint[1] = (1 - prob_joint[0] - prob_joint[3])/2;
    // P(1,0|S)
    prob_joint[2] = prob_joint[1];
    //Rcout<<prob_joint[0]<<"\t"<<prob_joint[1]<<"\t"<<prob_joint[2]<<"\t"<<prob_joint[3]<<std::endl;
    Pjoint_array.push_back(prob_joint);
    
    // add to loglikelihood 
    for(int pj_ind=0; pj_ind < 4; ++pj_ind){
      log_likelihood += emp_joint_counts[S_idx][pj_ind] * log(prob_joint[pj_ind]);
    }
    //std::cout<<"S="<<S<<"; loglik="<<log_likelihood<<std::endl;
  }
  
  
  
  
  //// GRADIENT CALCULATION  ///////
  //std::cout<<"Start calc FW_Jakobian"<<std::endl;
  std::vector<std::vector<double > > FW_Jakobian(max_spacing + 1,
                                                 std::vector<double >(n_ftp,0));
  std::vector<std::vector<double > > sigma_cumul_Jakobian(max_spacing + 1,
                                                          std::vector<double >(n_ftp,0));
  
  //   omp_set_num_threads(20);
  //   //std::cout<<"Number of threads="<<omp_get_max_threads();
  //   int ftp_nu=0;
  // #pragma omp parallel private(ftp_nu)
  // {
  //   
  // #pragma omp for schedule(dynamic)  
  for(int ftp_nu=0; ftp_nu < n_ftp;++ftp_nu){
    
    /// GRAD FOR PARTITION SUM AND SIGMA
    for(int d=1; d <= max_spacing; ++d){
      double fw_sum_tmp = 0;
      for(int ftp=0; ftp < n_ftp; ++ftp){
        if(d - ftp_lengths[ftp] < 0)
          break;
        fw_sum_tmp += FW[d - ftp_lengths[ftp]] * ftp_start_priors[ftp]/(R_const * ftp_lengths[ftp_nu]) - 
          FW_Jakobian[d - ftp_lengths[ftp]][ftp_nu] * ftp_start_priors[ftp];
      }
      if(d - ftp_lengths[ftp_nu] >=0){
        FW_Jakobian[d][ftp_nu] += FW[d - ftp_lengths[ftp_nu]]/(R_const * ftp_lengths[ftp_nu]);
      }
      FW_Jakobian[d][ftp_nu] -= fw_sum_tmp;
      
      sigma_cumul_Jakobian[d][ftp_nu] = sigma_cumul_Jakobian[d-1][ftp_nu] + 
        (FW[d-1]/R_const) * (kronecker_delta(0,ftp_nu)  - ftp_start_priors[0]/ftp_lengths[ftp_nu]) + 
        ftp_start_priors[0] * FW_Jakobian[d-1][ftp_nu];
      
    }
    
    // GRAD FOR LOGLIK FUNCTION//
    for(int S_idx = 0; S_idx < n_spac; ++S_idx){
      int S = spacings[S_idx];
      double sigma_sum0_grad = 0;
      double sigma_sum1_grad = 0;
      for(int ftp=0; ftp < n_ftp; ++ftp){
        int ftp_len = ftp_lengths[ftp];
        //std::cout<<"S="<<S<<"; ftp_nu="<<ftp_nu<<"; ftp="<<ftp<<"; ftp_len="<<ftp_len<<std::endl;
        if(S - ftp_len - 1 <=0)
          break;
        //std::cout<<"S="<<S<<"; ftp_nu="<<ftp_nu<<"; ftp="<<ftp<<"; sigma_cumul_Jakobian[S - ftp_len - 1]"<<sigma_cumul[S - ftp_len - 1]<<std::endl;
        double tmp_sum = ftp_cover_priors(ftp) * sigma_cumul_Jakobian[S - ftp_len - 1][ftp_nu]/ftp_len;
        sigma_sum0_grad += alpha0_vec[ftp] * tmp_sum;
        sigma_sum1_grad += alpha1_vec[ftp] * tmp_sum;
        
      }
      if(S - ftp_lengths[ftp_nu]  - 1 >= 0){
        sigma_sum0_grad += (alpha0_vec[ftp_nu]/ftp_lengths[ftp_nu]) * sigma_cumul[S - ftp_lengths[ftp_nu] - 1];
        sigma_sum1_grad += (alpha1_vec[ftp_nu]/ftp_lengths[ftp_nu]) * sigma_cumul[S - ftp_lengths[ftp_nu] - 1];
      }
      
      std::vector<double > prob_joint_grad(4,0);
      //// dP(0,0|S)/d(ftp_nu)
      prob_joint_grad[0] = alpha0 * (beta0 - alpha0) * kronecker_delta(0,ftp_nu) + (beta0 - alpha0) * (
        (R_const * alpha0 + (beta0 - alpha0) * ftp_cover_priors(0)) * sigma_cumul_Jakobian[S - 1][ftp_nu] + 
          sigma_cumul[S-1] * (alpha0/ftp_lengths[ftp_nu] + (beta0 - alpha0) * kronecker_delta(0,ftp_nu)) - 
          sigma_sum0_grad
      );
      
      
      // dP(0,1|S)
      prob_joint_grad[1] = alpha1 * (beta0 - alpha0) * kronecker_delta(0,ftp_nu) + (beta1 - alpha1) * (
        (R_const * alpha0 + (beta0 - alpha0) * ftp_cover_priors(0)) * sigma_cumul_Jakobian[S - 1][ftp_nu] + 
          sigma_cumul[S-1] * (alpha0/ftp_lengths[ftp_nu] + (beta0 - alpha0) * kronecker_delta(0,ftp_nu)) - 
          sigma_sum0_grad
      );  
      
      // dP(1,0|S)
      prob_joint_grad[2] = alpha0 * (beta1 - alpha1) * kronecker_delta(0,ftp_nu) + (beta0 - alpha0) * (
        (R_const * alpha1 + (beta1 - alpha1) * ftp_cover_priors(0)) * sigma_cumul_Jakobian[S - 1][ftp_nu] + 
          sigma_cumul[S-1] * (alpha1/ftp_lengths[ftp_nu] + (beta1 - alpha1) * kronecker_delta(0,ftp_nu)) -
          sigma_sum1_grad
      );
      
      
      // dP(1,1|S)
      prob_joint_grad[3] = alpha1 * (beta1 - alpha1) * kronecker_delta(0,ftp_nu) + (beta1 - alpha1) * (
        (R_const * alpha1 + (beta1 - alpha1) * ftp_cover_priors(0)) * sigma_cumul_Jakobian[S - 1][ftp_nu] + 
          sigma_cumul[S-1] * (alpha1/ftp_lengths[ftp_nu] + (beta1 - alpha1) * kronecker_delta(0,ftp_nu)) - 
          sigma_sum1_grad
      );
      
      
      
      // add to loglikelihood gradient
      
      //std::cout<<"log_likelihood_gradient.size()="<<log_likelihood_gradient.size()<<std::endl;
      for(int pj_ind=0; pj_ind < 4; ++pj_ind){
        log_likelihood_gradient[ftp_nu] += (emp_joint_counts[S_idx][pj_ind]/Pjoint_array[S_idx][pj_ind]) * prob_joint_grad[pj_ind];
      }
      
      
    }
    
    
  }
  //}  
  //std::cout<<"END of function"<<std::endl;
  
  // std::cout<<std::endl;
  // for(int ftp_nu=0;ftp_nu<n_ftp;++ftp_nu){
  //   std::cout<<log_likelihood_gradient[ftp_nu]<<",";
  // }
  // std::cout<<std::endl;
  
  return(log_likelihood);
  
  
}

}







#endif