#ifndef _bindobjmodel_h_
#define _bindobjmodel_h_

#include "DNAbinding_object-class.hpp"
//#include "nucleosomemodel.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include "parameters-class.hpp"
#include "utils_globvars.hpp"
#include "SMFdataset-class.hpp"
using namespace std;


class binding_object_model:public DNAbinding_object
{
 public:
  vector<vector<double> > mat;
  vector<vector<double> > normmat;
  vector<vector<double> > firstLastRatios; //  matrix containing ratios of first and last positions within WM
                                           // for all combintations of letters
                                           //                         letter at first WM position
                                           //                               | 0 | 1 | 2 |
                                           //  letter at last WM position 0 |...|...|...|
                                           //                             1 |...|...|...|
                                           //                             2 |...|...|...|
  binding_object_model(const vector<double > &_protect_prob,
                       const double _prior,
                       const string _name,
                       const string _group);

  virtual ~binding_object_model();
  virtual void print() const;
  virtual void print_normalized() const;

  virtual double get_score(const SMFdataset& SEQUENCES,
                           int seq,
                           int position) const;
  // method to pre-calculate footprint scores for a given molecule
  virtual vector<double > get_seq_scores_vec(const fragProtectData& fragData) const;
  void normalize();

};


#endif
