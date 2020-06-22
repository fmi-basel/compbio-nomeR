#ifndef _ftp_model_loglik_and_grad_templates_hpp_
#define _ftp_model_loglik_and_grad_templates_hpp_

#include "ftp_model_loglik_grad_v1.hpp"


// templates for functions to use in stan
// currently not used anywhere yet



//using stan::math::var;

template <typename T0__, typename T2__, typename T3__>
stan::promote_args_t<T0__, T2__,T3__>
ftp_model_lg_pr(const Eigen::Matrix<T0__, -1, 1>& ftp_cover_priors,
                const int& n_ftp, const T2__& bg_protect_prob,
                const T3__& footprint_protect_prob,
                const std::vector<int>& spacings, const int& n_spac,
                const int& max_spacing,
                const std::vector<std::vector<int>>& emp_joint_counts,
                std::ostream* pstream__) {
  throw std::logic_error("not implemented");  // this should never be called
}


//overloaded with double
//example from https://discourse.mc-stan.org/t/what-is-the-best-way-to-add-a-user-function-that-is-used-in-likelihood-calculation/4126/5
// template<>
// double spec_fun<double, double>(const double& p, const Eigen::Matrix<double,-1,1>& x, std::ostream* pstream__) {
//   return linear_sde::spec_fun(p, x);
// }

// template to <double,double,double>
template<>
double ftp_model_lg_pr<double, double, double>(const Eigen::Matrix<double, -1, 1>& ftp_cover_priors,  
                                               const int& n_ftp, // number of footprint. must equal to ftp_cover_priors.size()
                                               const double& bg_protect_prob, // emission prob for bg
                                               const double& footprint_protect_prob, // emission pros for footprints
                                               const std::vector<int >& spacings, // vector containing of spacing values
                                               const int& n_spac, // lengths of spacings vector. must be equal to spacings.size()
                                               const int& max_spacing, // maximum spacing in vector spacings. must be equal to max(spacings)
                                               const std::vector<std::vector<int > >& emp_joint_counts,
                                               // columns must be N00, N01, N10, N11
                                               std::ostream* pstream__
){
  // example from https://github.com/stan-dev/math/wiki/Adding-a-new-function-with-known-gradients
  // vector<var> foo(const vector<var>& x) {
  //   vector<double> a = value_of(x);
  //   double fa = bar::foo(a);
  //   vector<double> grad_fa = bar::grad_foo(a);  
  //   return precomputed_gradients(fa, x, grad_fa);
  // }
  
  // calculate lp and grad
  std::vector<double > log_likelihood_gradient(n_ftp,0);
  //Eigen::VectorXd log_likelihood_gradient = Eigen::VectorXd::Zero(n_ftp);
  //std::cout<<"In <doubel, double, double>"<<std::endl;
  double loglik = ftpmodfunctions::ftp_model_loglik_grad(ftp_cover_priors,
                                                         n_ftp,
                                                         bg_protect_prob,
                                                         footprint_protect_prob,
                                                         spacings, 
                                                         n_spac,
                                                         max_spacing,
                                                         emp_joint_counts, 
                                                         log_likelihood_gradient,
                                                         pstream__
                                                           
  );
  
  
  return(loglik);
  
}



// overloaded with var

// template for <var, double, double>
template<>
stan::math::var ftp_model_lg_pr<var, double, double>(const Eigen::Matrix<stan::math::var, -1, 1>& ftp_cover_priors,   
                                         const int& n_ftp, // number of footprint. must equal to ftp_cover_priors.size()
                                         const double& bg_protect_prob, // emission prob for bg
                                         const double& footprint_protect_prob, // emission pros for footprints
                                         const std::vector<int >& spacings, // vector containing of spacing values
                                         const int& n_spac, // lengths of spacings vector. must be equal to spacings.size()
                                         const int& max_spacing, // maximum spacing in vector spacings. must be equal to max(spacings)
                                         const std::vector<std::vector<int > >& emp_joint_counts,
                                         // columns must be N00, N01, N10, N11
                                         std::ostream* pstream__
){
  // example from https://github.com/stan-dev/math/wiki/Adding-a-new-function-with-known-gradients
  // vector<var> foo(const vector<var>& x) {
  //   vector<double> a = value_of(x);
  //   double fa = bar::foo(a);
  //   vector<double> grad_fa = bar::grad_foo(a);  
  //   return precomputed_gradients(fa, x, grad_fa);
  // }
  
  
  // extract vector of priors
  Eigen::Matrix<double, -1, 1> ftp_cover_priors_ = value_of(ftp_cover_priors);
  
  // calculate lp and grad
  std::vector<double > log_likelihood_gradient(n_ftp,0);
  double loglik = ftpmodfunctions::ftp_model_loglik_grad(ftp_cover_priors_,
                                                         n_ftp,
                                                         bg_protect_prob,
                                                         footprint_protect_prob,
                                                         spacings, 
                                                         n_spac,
                                                         max_spacing,
                                                         emp_joint_counts, 
                                                         log_likelihood_gradient,
                                                         pstream__
                                                           
  );
  
  return precomputed_gradients(loglik, ftp_cover_priors, log_likelihood_gradient);
  
  
}


#endif