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
#include "utils_globvars.hpp"
//#include "yaml-cpp/yaml.h"
using namespace std;


class parameters
{
    public:
      string printoutonly;// will be always "All" in R wrapper

      double bgcoverprob; // emission probability of a protected position in background (open/accessible) chromatin
      double bgprior; // prior probability of background
      parameters();
      ~parameters();


      // simple constructor with  all parameters
      parameters(double _bgcoverprob,
                 double _bgprior);
      parameters(const Rcpp::NumericVector& _bgprotectprob,
                 const Rcpp::NumericVector& _bgprior);

      void setParams(double _bgcoverprob,
                     double _bgprior);

      parameters & operator = (const parameters & other);
      void print();

      void clear();
};



#endif
