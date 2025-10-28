#ifndef _predict_hpp_
#define _predict_hpp_

#include "parameters-class.hpp"
#include "utils_globvars.hpp"
#include "DNAbindobj_vector-class.hpp"
#include "fragProtectData-class.hpp"
#include "SMFdataset-class.hpp"
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
	
	const parameters& PARAMS; // reference to an object containing parameters
	const DNAbind_obj_vector& BINDING_OBJECTS; // reference to an object containing vector of footprint models as well as background model
	const SMFdataset& SEQUENCES; // reference to an object containing SMF data

  // maps of footprint indices to data within Predict
  vector<int > print_indexes;	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
  vector<vector<int > > names2indexes; // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  vector<string > print_names; // this array contain names of the objects that will be printed
  vector<vector<int > > names2indicesinprobarray; // this array contains map i - index in names to subarray of indices in Prob array


public:
  // constructors/destructor
  //Predict();
	Predict(const SMFdataset& refSmfData,
         const DNAbind_obj_vector& refFtp_models,
         const parameters& refParams);

  bool Create();

  ~Predict();
  void clear();

  Rcpp::List calcStartCoverProbs(bool report_prediction_in_flanks,
                                 int ncpu);

};



#endif
