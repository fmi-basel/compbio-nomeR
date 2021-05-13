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
#include "refSeqInfo-class.hpp"
#include "coocCntTable-class.hpp"
#include "regionData-class.hpp"
#include "protectStats-class.hpp"
#include "predict-class.hpp"
#include "fetch_data_from_bam.hpp"
#include "nomeseqdata.h"

using namespace std;


bool _VERBOSE_ = 0;



// [[Rcpp::export]]
Rcpp::List fetch_cooc_ctable_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::IntegerVector& max_spac,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax);


// [[Rcpp::export]]
Rcpp::List fetch_data_matrix_from_bams_cpp(const Rcpp::CharacterVector& whichContext,
                                           const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax);

// [[Rcpp::export]]
Rcpp::List fetch_protect_stats_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                             const Rcpp::CharacterVector& regionChr,
                                             const Rcpp::IntegerVector& regionStart,
                                             const Rcpp::IntegerVector& regionEnd,
                                             const Rcpp::CharacterVector& seqstring,
                                             const Rcpp::IntegerVector& seqStart,
                                             const Rcpp::IntegerVector& seqEnd,
                                             const Rcpp::LogicalVector& remove_nonunique,
                                             const Rcpp::IntegerVector& clip_until_nbg,
                                             const Rcpp::NumericVector& max_protect_frac,
                                             const Rcpp::NumericVector& max_bisC_meth,
                                             const Rcpp::IntegerVector& min_bisC_size,
                                             const Rcpp::IntegerVector& mapqMin,
                                             const Rcpp::IntegerVector& mapqMax);


// [[Rcpp::export]]
Rcpp::List run_cpp_nomeR(const Rcpp::List& data,
                         const Rcpp::CharacterVector& fragnames,
                         const Rcpp::List& binding_models,
                         const Rcpp::NumericVector& bgprotectprob,
                         const Rcpp::NumericVector& bgprior,
                         const Rcpp::LogicalVector& report_prediction_in_flanks,
                         const Rcpp::NumericVector& Ncpu,
                         const Rcpp::LogicalVector& verbose);

// [[Rcpp::export]]
Rcpp::List count_spacing_freq_cpp(const Rcpp::List& data,
                                  const Rcpp::CharacterVector& fragnames,
                                  const Rcpp::IntegerVector& maxspacing,
                                  const Rcpp::IntegerVector& maxwmlen);


// [[Rcpp::export]]
Rcpp::List calculate_theor_joint_prob_cpp(const Rcpp::NumericVector& ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                          const Rcpp::NumericVector& bg_protect_prob,
                                          const Rcpp::NumericVector& footprint_protect_prob,
                                          const Rcpp::IntegerVector& max_spacing);



#endif