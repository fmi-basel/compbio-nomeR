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
                      const IntegerVector& priorEM_max_steps
                      ) {

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
  
  PARAMS.setParams(bgcoverprob_,
                   bgprior_,
                   bound_fit_tol_,
                   bound_min_fderiv_val_,
                   bound_max_steps_,
                   priorEM_fit_tol_,
                   priorEM_max_steps_);
  
  //PARAMS.print();
  
  // create object with binding models
  //cout<<"Size "<<binding_models.size()<<endl;
  extern DNAbind_obj_vector BINDING_OBJECTS;
  BINDING_OBJECTS.create(binding_models,PARAMS);
  //BINDING_OBJECTS.print();
  //cout<<"MAX WM length "<<BINDING_OBJECTS.maxwmlen<<endl;
  // create object with NOMeSeq data
  extern NOMeSeqData SEQUENCES;
  SEQUENCES.create(data,
                   BINDING_OBJECTS.maxwmlen);
  //SEQUENCES.PrintNames();
  
  extern Forward_Backward_algorithm FORWARD_BACKWARD;
  FORWARD_BACKWARD.Create();
  if(run_priorEM_){
    FORWARD_BACKWARD.Run_priorEM();
  } else {
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
                        int maxspacing,
                        int maxwmlen = 0){
  NOMeSeqData nome_data;
  nome_data.create(data,
         maxwmlen);
  
  List spac_freq = nome_data.R_export_spacing_freq(maxspacing);
  return spac_freq;
}
