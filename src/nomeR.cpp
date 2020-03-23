#include "parameters.h"
#include "nomeseqdata.h"
#include "DNAbindobj_vector.h"
#include "forward_backward.h"
#include <Rcpp.h>

using namespace Rcpp;


// global objects
parameters PARAMS;
DNAbind_obj_vector BINDING_OBJECTS;
NOMeSeqData SEQUENCES;
Forward_Backward_algorithm FORWARD_BACKWARD;
int _NCPU_ = 1;


// [[Rcpp::export]]
List run_cpp_nomeR(const List& data,
                      const List& binding_models,
                      const NumericVector& bgprotectprob,
                      const NumericVector& bgprior,
                      const NumericVector& bound_fit_tol,
                      const NumericVector& bound_min_fderiv_val,
                      const IntegerVector& bound_max_steps,
                      const LogicalVector& run_priorEM,
                      const NumericVector& priorEM_fit_tol,
                      const IntegerVector& priorEM_max_steps,
                      const NumericVector& Ncpu,
                      const LogicalVector& verbose
                      ) {
  // set number of threads
  
  
  extern int _NCPU_;
  _NCPU_ = as<int >(Ncpu);
  if(verbose){
    Rcout<<"Input Ncpu "<<Ncpu<<endl;
    Rcout<<"Set _NCPU_ "<<_NCPU_<<endl;
  }
  
  
  
  // create object with parameters
  extern parameters PARAMS;
  double bgcoverprob_ = as<double >(bgprotectprob);
  double bgprior_ = as<double >(bgprior);
  double bound_fit_tol_ = as<double >(bound_fit_tol);
  double bound_min_fderiv_val_ = as<double >(bound_min_fderiv_val);
  int bound_max_steps_ = as<int >(bound_max_steps);
  
  bool run_priorEM_ = as<bool >(run_priorEM);
  double priorEM_fit_tol_ = as<double >(priorEM_fit_tol);
  int priorEM_max_steps_ = as<int >(priorEM_max_steps);
  
  if(verbose){
    Rcout<<"Creating PARAMS object..."<<endl;
  }
  PARAMS.setParams(bgcoverprob_,
                   bgprior_,
                   bound_fit_tol_,
                   bound_min_fderiv_val_,
                   bound_max_steps_,
                   priorEM_fit_tol_,
                   priorEM_max_steps_);
  
  //PARAMS.print();
  
  // create object with binding models
  
  
  if(verbose){
    Rcout<<"Creating BINDING_OBJECTS object..."<<endl;
  }
  extern DNAbind_obj_vector BINDING_OBJECTS;
  BINDING_OBJECTS.create(binding_models,PARAMS);
  //BINDING_OBJECTS.print();
  //Rcout<<"MAX WM length "<<BINDING_OBJECTS.maxwmlen<<endl;
  // create object with NOMeSeq data
  extern NOMeSeqData SEQUENCES;
  
  if(verbose){
    Rcout<<"Creating SEQUENCES object..."<<endl;
  }
  SEQUENCES.create(data,
                   BINDING_OBJECTS.maxwmlen);
  //SEQUENCES.PrintNames();
  
  extern Forward_Backward_algorithm FORWARD_BACKWARD;
  if(verbose){
    Rcout<<"Creating FORWARD_BACKWARD object..."<<endl;
  }
  FORWARD_BACKWARD.Create();
  
  
  if(run_priorEM_){
    if(verbose){
      Rcout<<"Running FORWARD_BACKWARD.Run_priorEM()..."<<endl;
    }
    FORWARD_BACKWARD.Run_priorEM();
  } else {
    if(verbose){
      Rcout<<"Running FORWARD_BACKWARD.Run()..."<<endl;
    }
    FORWARD_BACKWARD.Run();
  }
  
  List startProbdf = FORWARD_BACKWARD.getStartProbDF();
  List coverProb = FORWARD_BACKWARD.getCoverProbDF();
  List genSummary = FORWARD_BACKWARD.getGenomeSummaryDF();
  List output_data = List::create( Named("START_PROB") = startProbdf,
                         Named("COVER_PROB") = coverProb,
                         Named("SUMMARY") = genSummary);
  
  PARAMS.clear();
  BINDING_OBJECTS.clear();
  SEQUENCES.clear();
  FORWARD_BACKWARD.clear();
  return output_data;
}

// [[Rcpp::export]]
List count_spacing_freq_cpp(const List& data,
                            const IntegerVector& maxspacing,
                            const IntegerVector& maxwmlen){
  
  Rcout<<"In count_spacing_freq_cpp"<<endl;
  Rcout<<"Creting object for nome data 0"<<endl;
  NOMeSeqData nome_data;
  Rcout<<"Creting object for nome data 1"<<endl;
  int maxwmlen_ = as<int >(maxwmlen);
  int maxspacing_ = as<int >(maxspacing);
  nome_data.create(data,
         maxwmlen_);
  Rcout<<"Counting spacings"<<endl;
  List spac_freq = nome_data.R_export_spacing_freq(maxspacing_);
  // return spac_freq;
  Rcpp::List out_list;
  return out_list;
  
}
