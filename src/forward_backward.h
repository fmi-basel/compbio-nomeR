#ifndef _forward_backward_h_
#define _forward_backward_h_

#include "parameters.h"
#include "usefullfunctions.h"
#include "DNAbindobj_vector.h"
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

class Forward_Backward_algorithm
{

 public:
  vector<vector<double> > F;
  vector<vector<double> > R;
  vector<vector<vector<double> > > Prob;
  vector<vector<double > > genomesummary; //this vector contains expected prior, expected number of sites and expected coverage for the whole genome and for each TF
  
  vector<int > print_indexes;	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
  vector<vector<int > > names2indexes; // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  vector<string > print_names; // this array contain names of the objects that will be printed

  vector<vector<int > > names2indicesinprobarray; // this array contains map i - index in names to subarray of indices in Prob array



  Forward_Backward_algorithm();
  bool Create();
  ~Forward_Backward_algorithm();
  
  void Run(vector<double > priors,
           int ncpu);
  void Run(int ncpu);

  vector<double > getPriors();
  
  
  
  Rcpp::List getStartProbDF();
  Rcpp::List getCoverProbDF();
  Rcpp::List getGenomeSummaryDF();
  

  void SetGenomeSummary();
  double get_coverage_prob_at_pos_name_index(int seq, int pos, int name_index);

  double get_coverage_prob_at_pos(int seq, int pos, int wm); // wm here is an index in Prob
  
  void clear();
};


#endif
