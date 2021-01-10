#ifndef _predict_footprints_h_
#define _predict_footprints_h_

#include "parameters.h"
#include "usefullfunctions.h"
#include "DNAbindobj_vector.h"
#include "nomeseqdata.h"
#include "forward_backward.h"
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

class predict_footprints
{
  parameters PARAMS;
  DNAbind_obj_vector BINDING_OBJECTS;
  NOMeSeqData SEQUENCES;
  Forward_Backward_algorithm FORWARD_BACKWARD;
  //int _NCPU_ = 1;
  bool _VERBOSE_ = 0;
  
public:
  Forward_Backward_algorithm();
  ~Forward_Backward_algorithm();
  void Run(int ncpu);
};

#endif
