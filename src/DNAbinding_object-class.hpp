/* This abstract class defines any objects wich can bind to DNA
   Weight matrix and nucleosome will be descendants of this class */

#ifndef _DNAbinding_object_hpp_
#define _DNAbinding_object_hpp_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "SMFdataset-class.hpp"
using namespace std;

class DNAbinding_object{
 public:
  string classname;
  string name;
  string group;
  int len;
  double prior;
  double initialprior;

  DNAbinding_object();
  virtual ~DNAbinding_object() = 0;

  virtual double get_score(const SMFdataset& SEQUENCES,
                           int seq,
                           int position) const = 0;
  // virtual method to pre-calculate footprint scores for a given molecule
  virtual vector<double > get_seq_scores_vec(const fragProtectData& fragData) const = 0;
  virtual void print() const = 0;
  virtual void print_normalized() const=0;
};

#endif
