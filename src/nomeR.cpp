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
bool _VERBOSE_ = 0;

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
  
  //set verbose
  extern bool _VERBOSE_;
  _VERBOSE_ = as<bool >(verbose);
  
  // set number of threads
  extern int _NCPU_;
  _NCPU_ = as<int >(Ncpu);
  if(_VERBOSE_){
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
  
  if(_VERBOSE_){
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
  
  
  if(_VERBOSE_){
    Rcout<<"Creating BINDING_OBJECTS object..."<<endl;
  }
  extern DNAbind_obj_vector BINDING_OBJECTS;
  BINDING_OBJECTS.create(binding_models,PARAMS);
  //BINDING_OBJECTS.print();
  //Rcout<<"MAX WM length "<<BINDING_OBJECTS.maxwmlen<<endl;
  // create object with NOMeSeq data
  extern NOMeSeqData SEQUENCES;
  
  if(_VERBOSE_){
    Rcout<<"Creating SEQUENCES object..."<<endl;
  }
  SEQUENCES.create(data,
                   BINDING_OBJECTS.maxwmlen);
  //SEQUENCES.PrintNames();
  
  extern Forward_Backward_algorithm FORWARD_BACKWARD;
  if(_VERBOSE_){
    Rcout<<"Creating FORWARD_BACKWARD object..."<<endl;
  }
  FORWARD_BACKWARD.Create();
  
  
  if(run_priorEM_){
    if(_VERBOSE_){
      Rcout<<"Running FORWARD_BACKWARD.Run_priorEM()..."<<endl;
    }
    FORWARD_BACKWARD.Run_priorEM();
  } else {
    if(_VERBOSE_){
      Rcout<<"Running FORWARD_BACKWARD.Run()..."<<endl;
    }
    FORWARD_BACKWARD.Run();
  }
  if(_VERBOSE_){
    Rcout<<"Running FORWARD_BACKWARD.getStartProbDF()..."<<endl;
  }
  List startProbdf = FORWARD_BACKWARD.getStartProbDF();
  if(_VERBOSE_){
    Rcout<<"Running FORWARD_BACKWARD.getCoverProbDF()..."<<endl;
  }
  List coverProb = FORWARD_BACKWARD.getCoverProbDF();
  if(_VERBOSE_){
    Rcout<<"Running FORWARD_BACKWARD.getGenomeSummaryDF()..."<<endl;
  }
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
  
  int maxwmlen_ = as<int >(maxwmlen);
  int maxspacing_ = as<int >(maxspacing);
  
  NOMeSeqData nome_data;
  nome_data.create(data,maxwmlen_);
  vector<vector<int> > freq_mat = nome_data.count_freq_for_spacings(maxspacing_);
  vector<int > spacings;
  vector<int > freq00;
  vector<int > freq01;
  vector<int > freq02;
  vector<int > freq10;
  vector<int > freq11;
  vector<int > freq12;
  
  vector<int > freq20;
  vector<int > freq21;
  vector<int > freq22;
  for(int i=0; i<freq_mat.size(); ++i){
    spacings.push_back(freq_mat[i][0]);
    freq00.push_back(freq_mat[i][1]);
    freq01.push_back(freq_mat[i][2]);
    freq02.push_back(freq_mat[i][3]);
    freq10.push_back(freq_mat[i][4]);
    freq11.push_back(freq_mat[i][5]);
    freq12.push_back(freq_mat[i][6]);
    freq20.push_back(freq_mat[i][7]);
    freq21.push_back(freq_mat[i][8]);
    freq22.push_back(freq_mat[i][9]);
    
  }
  
  List export_list;
  
  export_list.push_back(Rcpp::wrap(spacings),"S");
  export_list.push_back(Rcpp::wrap(freq00),"N00");
  export_list.push_back(Rcpp::wrap(freq01),"N01");
  export_list.push_back(Rcpp::wrap(freq10),"N10");
  export_list.push_back(Rcpp::wrap(freq11),"N11");
  export_list.push_back(Rcpp::wrap(freq02),"N0na");
  export_list.push_back(Rcpp::wrap(freq12),"N1na");
  export_list.push_back(Rcpp::wrap(freq20),"Nna0");
  export_list.push_back(Rcpp::wrap(freq21),"Nna1");
  export_list.push_back(Rcpp::wrap(freq22),"Nnana");
  
  
  nome_data.clear();
  return export_list;
  
}

// [[Rcpp::export]]
List calculate_theor_joint_prob_cpp(const NumericVector& ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                const NumericVector& bg_protect_prob,
                                const NumericVector& footprint_protect_prob,
                                const IntegerVector& max_spacing){
  
  vector<double > ftp_cover_priors_ = Rcpp::as<vector<double > >(ftp_cover_priors);
  double bg_protect_prob_ = as<double >(bg_protect_prob);
  double footprint_protect_prob_ = as<double >(footprint_protect_prob);
  int max_spacing_ = as<int >(max_spacing);
  
  
  DNAbind_obj_vector ftp_array;
  vector<vector<double > > p_joint = ftp_array.calc_theor_joint_prob(ftp_cover_priors_,
                                                                     bg_protect_prob_,
                                                                     footprint_protect_prob_,
                                                                     max_spacing_);
  
  // create Rcpp object
  
  vector<int > spacings;
  vector<double > p00;
  vector<double > p01;
  vector<double > p10;
  vector<double > p11;
  
  
  for(int i = 0; i < p_joint.size(); ++i){
    //Rcpp::Rcout<<"S="<<p_joint[i][0]<<"; "<<p_joint[i][1]<<"; "<<p_joint[i][2]<<"; "<<p_joint[i][3]<<"; "<<p_joint[i][4]<<endl;
    spacings.push_back(p_joint[i][0]);
    p00.push_back(p_joint[i][1]);
    p01.push_back(p_joint[i][2]);
    p10.push_back(p_joint[i][3]);
    p11.push_back(p_joint[i][4]);
  }
  
  Rcpp::List export_list;
  export_list.push_back(Rcpp::wrap(spacings),"S");
  export_list.push_back(Rcpp::wrap(p00),"P00");
  export_list.push_back(Rcpp::wrap(p01),"P01");
  export_list.push_back(Rcpp::wrap(p10),"P10");
  export_list.push_back(Rcpp::wrap(p11),"P11");
  
  
  ftp_array.clear();
  return export_list;
}

