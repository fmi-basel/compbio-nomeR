#ifndef _dnabindobjvec_h_
#define _dnabindobjvec_h_
#include "DNAbinding_object-class.hpp"
#include <Rcpp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include "parameters-class.hpp"
#include <fstream>
#include "binding_object_model-class.hpp"
#include "background-class.hpp"
#include "utils_globvars.hpp"

using namespace std;


class DNAbind_obj_vector{
 vector<DNAbinding_object *> objvector;
 int size;
 public:
   int numberofobjects;
   int maxwmlen;
   vector<string > names;
   vector<vector<int > > names2index;
   vector<double > priors;
   
   
   
  int Size();
  DNAbind_obj_vector();
  ~DNAbind_obj_vector();
  //int create(parameters &params,NOMeSeqData &sequences);
  
  int create(Rcpp::List _bind_objs,
             parameters &params);
  
  DNAbinding_object* operator [](int i);
  
  void print();
  
  void clear();
  
  vector<vector<double > > calc_theor_joint_prob(vector<double > ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector
                                                                                    // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
                                                 double bg_protect_prob,
                                                 double footprint_protect_prob,
                                                 int max_spacing);
  // Rcpp::List R_export_calc_theor_joint_prob(vector<double > ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector
  //                                           // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
  //                                           double bg_protect_prob,
  //                                           double footprint_protect_prob,
  //                                           int max_spacing);
  
};






#endif
