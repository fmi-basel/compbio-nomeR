#ifndef _background_h_
#define _background_h_

#include "DNAbinding_object-class.hpp"
//#include "nucleosomemodel.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include "utils_globvars.hpp"
#include "parameters-class.hpp"
#include "nomeseqdata.h"
using namespace std;

class Background:public DNAbinding_object
{
 public:
  double bgcoverprob;

  vector<double>  bgmodel;
  Background(parameters &params);
  virtual ~Background();
  virtual void print() const;
  
  virtual double get_score(NOMeSeqData& SEQUENCES,
                           int seq,
                           int position)const;
  virtual void print_normalized() const;
};

#endif
