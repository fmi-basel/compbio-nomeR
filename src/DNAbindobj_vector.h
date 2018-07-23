#ifndef _dnabindobjvec_h_
#define _dnabindobjvec_h_
#include "DNAbinding_object.h"
#include <Rcpp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include "parameters.h"
#include <fstream>
#include "binding_object_model.h"
#include "background.h"
#include "usefullfunctions.h"

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
};

#endif
