#ifndef _bam_utils_export_rcpp_hpp_
#define _bam_utils_export_rcpp_hpp_

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
#include "fetch_data_from_bam.hpp"
#include "protectStats-class.hpp"



using namespace std;

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

// TODO

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
#endif