#include "predict.hpp"
#include <Rcpp.h>

using namespace Rcpp;

bool _VERBOSE_ = 0;

// [[Rcpp::export]]
List run_cpp_nomeR(const List& data,
                   const CharacterVector& fragnames,
                   const List& binding_models,
                   const NumericVector& bgprotectprob,
                   const NumericVector& bgprior,
                   const NumericVector& Ncpu,
                   const LogicalVector& verbose
) {
  
  //set verbose
  extern bool _VERBOSE_;
  _VERBOSE_ = as<bool >(verbose);
  
  int Ncpu_ = as<int >(Ncpu);
#ifndef _OPENMP
  Rcout<<"nomeR was compiled without OpenMP. ncpu does not have effect.\n";
#endif
  
  
  
  if(_VERBOSE_){
    Rcout<<"Creating Predict object..."<<endl;
  }
  Predict predict(data,
                  fragnames,
                  binding_models,
                  bgprotectprob,
                  bgprior);
  
  // run prediction
  if(_VERBOSE_){
    Rcout<<"Calculating posterior binding probabilities..."<<endl;
  }
  
  List output_data;
  if(predict.Run(Ncpu_)){
    //predict.Run(Ncpu_);
    
    if(_VERBOSE_){
      Rcout<<"Running predict.getStartProbDF()..."<<endl;
    }
    List startProbdf = predict.getStartProbDF();
    if(_VERBOSE_){
      Rcout<<"Running predict.getCoverProbDF()..."<<endl;
    }
    List coverProb = predict.getCoverProbDF();
    if(_VERBOSE_){
      Rcout<<"Running predict.getGenomeSummaryDF()..."<<endl;
    }
    List genSummary = predict.getGenomeSummaryDF();
    output_data = List::create( Named("START_PROB") = startProbdf,
                                Named("COVER_PROB") = coverProb,
                                Named("SUMMARY") = genSummary);
  } else {
    output_data = List::create( Named("START_PROB") = R_NilValue,
                                Named("COVER_PROB") = R_NilValue,
                                Named("SUMMARY") = R_NilValue);
    
    }
  
  _VERBOSE_ = 0;
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

