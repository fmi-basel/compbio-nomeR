#include "binding_object_model.h"

binding_object_model::binding_object_model(const vector<double > &_protect_prob,
                     const double _prior,
                     const string _name){
  
  classname = "wm";
  
  
  prior = _prior;
  initialprior = _prior;
  name = _name;
  
  
  for(int pos=0;pos < _protect_prob.size();++pos){
    vector<double > tmp;
    tmp.push_back(1-_protect_prob[pos]); // prob of unprotected
    tmp.push_back(_protect_prob[pos]);  // prob of protected
    tmp.push_back(1); //add 1 at the end to treat NAs in NOMe data  
    mat.push_back(tmp);
  }
  len=mat.size();
  normalize();
}




binding_object_model::~binding_object_model() {
 
}

void binding_object_model::print() const{
  Rcpp::Rcout <<"//\n";
  Rcpp::Rcout <<"NA\t"<< name <<"\n";

  Rcpp::Rcout <<"PRIOR\t"<< prior <<"\n";
  Rcpp::Rcout <<"POS\tPROT\tUNPROT\tNAs\n";
  for(int i=0;i<mat.size();i++){
    Rcpp::Rcout <<i+1;
    for(int j=0;j<mat[i].size();j++)
      Rcpp::Rcout << "\t"<< mat[i][j];
    Rcpp::Rcout <<"\n";
  }
  Rcpp::Rcout <<"//\n";
}


void binding_object_model::print_normalized() const{
  Rcpp::Rcout <<"//\n";
  Rcpp::Rcout <<"NA\t"<< name <<"\n";
  //Rcpp::Rcout <<"Orientation\t"<<orientation<<endl;
  Rcpp::Rcout <<"PRIOR\t"<< prior <<"\n";
  Rcpp::Rcout <<"POS\tPROT\tUNPROT\tNAs\n";
  for(int i=0;i<normmat.size();i++){
    Rcpp::Rcout <<i+1;
    for(int j=0;j<normmat[i].size();j++)
      Rcpp::Rcout << "\t"<< normmat[i][j];
    Rcpp::Rcout <<"\n";
  }
  Rcpp::Rcout <<"//\n";
}


double binding_object_model::get_score(int seq, int position) const{
  
  extern NOMeSeqData SEQUENCES;
  if(seq<0 || seq>=SEQUENCES.Size()){
    Rcpp::stop("binding_object_model::get_score: Index of sequence is out of range: ");
	// cerr<<"Wm::get_score: Index of sequence is out of range: "<<seq<<endl;
	// exit(1);
  }

  double score = 1;
  for(int i = position;i < position + len;++i){
    if(i >= 0 && i < SEQUENCES[seq].Size())
      score *= normmat[i - position][SEQUENCES[seq][i]];
  }
  return prior * score;
}


void binding_object_model::normalize(){
  double pseudocount = 0.0;
  vector<double > tmp(mat[0].size(),1);
  normmat.resize(mat.size(),tmp);
  for(int i=0;i<mat.size();++i){
    double summ=0;
    for(int j=0;j<=1;++j){
      summ+=mat[i][j] + pseudocount;
    }
    for(int j=0;j<=1;++j){
	    normmat[i][j] = (mat[i][j] + pseudocount)/summ;
    }
  }
}
