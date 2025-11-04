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
#include "SMFdataset-class.hpp"
using namespace std;

class Background:public DNAbinding_object
{
 public:
  double bgcoverprob;

  vector<double>  bgmodel;
  Background(const parameters& params);
  virtual ~Background();
  virtual void print() const;

  virtual double get_score(const SMFdataset& SEQUENCES,
                           int seq,
                           int position) const;
  // method to pre-calculate footprint scores for a given molecule
  virtual vector<double > get_seq_scores_vec(const fragProtectData& fragData) const;
  virtual void print_normalized() const;
};

#endif
