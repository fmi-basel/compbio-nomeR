#ifndef _functions_export_rcpp_hpp_
#define _functions_export_rcpp_hpp_


#include <Rcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cctype>
#include <vector>
#include <map>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdint>
#include <stdbool.h>

#include "parameters-class.hpp"
#include "SMFdataset-class.hpp"
#include "DNAbindobj_vector-class.hpp"
#include "predict-class.hpp"
using namespace std;


bool _VERBOSE_ = 0;



// [[Rcpp::export]]
Rcpp::List calcStartCoverProbs_cpp(const Rcpp::IntegerVector& fragIDs,     // vector with unique fragment IDs
                                   const Rcpp::IntegerVector& fragPos,     // vector with positions within each fragment, 1 - based!
                                   const Rcpp::IntegerVector& protectVec,  // vector with protection data, 0 - accessible; 1 - protected
                                   const Rcpp::List& binding_models,
                                   const Rcpp::NumericVector& bgprotectprob,
                                   const Rcpp::NumericVector& bgprior,
                                   const Rcpp::LogicalVector& report_prediction_in_flanks,
                                   const Rcpp::NumericVector& Ncpu,
                                   const Rcpp::LogicalVector& verbose);


// [[Rcpp::export]]
Rcpp::NumericMatrix count_spacing_freq_cpp(const Rcpp::IntegerVector& fragIDs,     // vector with unique fragment IDs
                                  const Rcpp::IntegerVector& fragPos,     // vector with positions within each fragment, 1 - based!
                                  const Rcpp::IntegerVector& protectVec,  // vector with protection data, 0 - accessible; 1 - protected
                                  const Rcpp::IntegerVector& maxspacing,
                                  const Rcpp::NumericVector& Ncpu,
                                  const Rcpp::LogicalVector& verbose);


// [[Rcpp::export]]
Rcpp::List calculate_theor_joint_prob_cpp(const Rcpp::NumericVector& ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                          const Rcpp::NumericVector& bg_protect_prob,
                                          const Rcpp::NumericVector& footprint_protect_prob,
                                          const Rcpp::IntegerVector& max_spacing);



#endif
