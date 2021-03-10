#ifndef _bindobjmodel_h_
#define _bindobjmodel_h_

#include "DNAbinding_object.h"
//#include "nucleosomemodel.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include "parameters.h"
#include "usefullfunctions.h"
#include "nomeseqdata.h"
using namespace std;


class binding_object_model:public DNAbinding_object
{
 public:
  vector<vector<double> > mat;
  vector<vector<double> > normmat;
  binding_object_model(const vector<double > &_protect_prob,
                       const double _prior,
                       const string _name);
  
  virtual ~binding_object_model();
  virtual void print() const;
  virtual void print_normalized() const;
  
  virtual double get_score(NOMeSeqData& SEQUENCES,
                           int seq,
                           int position)const;
  void normalize();

};


#endif
