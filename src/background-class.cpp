#include "background-class.hpp"

Background::Background(parameters &params){
  classname = "background";
  name = "background";
  len = 1;
  prior = params.bgprior;
  initialprior = prior;

  bgcoverprob = params.bgcoverprob;

  bgmodel.push_back(1-bgcoverprob);
  bgmodel.push_back(bgcoverprob);
  bgmodel.push_back(1); // this is for NAs
}


Background::~Background(){

}

void Background::print() const{
  Rcpp::Rcout << "Background parameters:\n";
  Rcpp::Rcout << "Background cover probability = "<<bgcoverprob<<endl;
  Rcpp::Rcout << "Background prior = "<<prior<<endl;
  
}
void Background::print_normalized() const{
  print();

}


double Background::get_score(NOMeSeqData& SEQUENCES,
                             int seq,
                             int position) const{
   
  //extern NOMeSeqData SEQUENCES;
  if(seq<0 || seq>=SEQUENCES.Size()){
    
  // Rcpp::Rcerr<<"Wm::get_score: Index of sequence is out of range: "<<seq<<endl;
    Rcpp::stop("Background::get_score: Index of sequence is out of range:");
  // exit(1);
  }

  double score = 1;
  for(int i = position;i < position + len;++i){
    if(i >= 0 && i < SEQUENCES[seq].Size())
      score *= bgmodel[SEQUENCES[seq][i]];
  }

  return prior * score;
  
}
