#ifndef _predict_hpp_
#define _predict_hpp_

#include "parameters-class.hpp"
#include "utils_globvars.hpp"
#include "DNAbindobj_vector-class.hpp"
#include "nomeseqdata.h"
#include "sequence.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <time.h>
#include <Rcpp.h>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::plugins(openmp)]]

class Predict
{
  vector<vector<double> > F; // forward partition sums for fragments
  vector<vector<double> > R; // backward partition sums for fragments
  vector<vector<vector<double> > > Prob; // start probabilities for footprints, fragments
  vector<vector<double > > genomesummary; //this vector contains expected prior, expected number of sites and expected coverage for the whole genome and for each TF
  
  vector<int > print_indexes;	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
  vector<vector<int > > names2indexes; // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  vector<string > print_names; // this array contain names of the objects that will be printed
  
  vector<vector<int > > names2indicesinprobarray; // this array contains map i - index in names to subarray of indices in Prob array
  
  parameters PARAMS; // object containing parameters
  DNAbind_obj_vector BINDING_OBJECTS; // object containing vector of footprint models as well as background model
  NOMeSeqData SEQUENCES; // NOMe-seq data sequences
  
public:
  // constructors/destructor
  Predict();
  Predict(const Rcpp::List& data,
          const Rcpp::CharacterVector& fragnames,
          const Rcpp::List& binding_models,
          const Rcpp::NumericVector& bgprotectprob,
          const Rcpp::NumericVector& bgprior);
  
  bool Create(const Rcpp::List& data,
              const Rcpp::CharacterVector& fragnames,
              const Rcpp::List& binding_models,
              const Rcpp::NumericVector& bgprotectprob,
              const Rcpp::NumericVector& bgprior);
  
  ~Predict();
  void clear();
  
  // run prediction
  bool Run(int ncpu);
  
  // utility functions
  vector<double > getPriors();
  void SetGenomeSummary();
  double get_coverage_prob_at_pos_name_index(int seq, int pos, int name_index);
  double get_coverage_prob_at_pos(int seq, int pos, int wm); // wm here is an index in Prob
  
  
  // output results
  Rcpp::List getStartProbDF(bool report_prediction_in_flanks);
  Rcpp::List getCoverProbDF();
  Rcpp::List getGenomeSummaryDF();
  
};



#endif