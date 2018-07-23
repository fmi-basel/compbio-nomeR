/* This abstract class defines any objects wich can bind to DNA 
   Weight matrix and nucleosome will be descendants of this class */

#ifndef _DNAbinding_object_h_
#define _DNAbinding_object_h_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
using namespace std;

class DNAbinding_object{
 public:
  string classname;
  string name;
  int len;
  double prior;
  double initialprior;

  DNAbinding_object();
  virtual ~DNAbinding_object() = 0;
  
  virtual double get_score(int seq, int position) const = 0;
  virtual void print() const = 0;
  virtual void print_normalized() const=0;
};

#endif
