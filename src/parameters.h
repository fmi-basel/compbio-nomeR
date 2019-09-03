#ifndef _parameters_h_
#define _parameters_h_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <Rcpp.h>
#include "usefullfunctions.h"
//#include "yaml-cpp/yaml.h"
using namespace std;


class parameters
{
    public:
      string printoutonly;// will be always "All" in R wrapper
      double BOUND_FIT_TOLERANCE;
      double BOUND_MIN_FDERIV_VAL;
      int BOUND_MAX_STEPS;

      double PRIOREM_FIT_TOLERANCE;
      int PRIOREM_MAX_STEPS;

      double bgcoverprob;
      double bgprior;
      parameters();
      ~parameters();
      
      
      // simple constructor with  all parameters
      parameters(double _bgcoverprob,
                 double _bgprior,
                 double _bound_fit_tol,
                 double _bound_min_fderiv_val,
                 int _bound_max_steps,
                 double _priorEM_fit_tol,
                 int _priorEM_max_steps);
      
      void setParams(double _bgcoverprob,
                     double _bgprior,
                     double _bound_fit_tol,
                     double _bound_min_fderiv_val,
                     int _bound_max_steps,
                     double _priorEM_fit_tol,
                     int _priorEM_max_steps);

      parameters & operator = (const parameters & other);
      void print();

      void clear();
};



#endif