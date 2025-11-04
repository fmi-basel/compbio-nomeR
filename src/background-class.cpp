#include "background-class.hpp"

Background::Background(const parameters& params){
  classname = "background";
  name = "background";
  group = "background";
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


double Background::get_score(const SMFdataset& SEQUENCES,
                             int seq,
                             int position) const{


  if(seq<0 || seq>=SEQUENCES.Size()){
    Rcpp::stop("Background::get_score: Index of sequence is out of range:");
  
  }

  double score = 1;
  for(int i = position; i < position + len;++i){
    if(i >= 0 && i < SEQUENCES[seq].Size())
      score *= bgmodel[SEQUENCES[seq][i]];
  }
  return prior * score;
}

vector<double > Background::get_seq_scores_vec(const fragProtectData& fragData) const{
	vector<double > scoresVec(fragData.Size(), prior);
	// as length of background is 1 we just fill the vector with prior * bgmodel[fragData[pos]]
	for(int pos = 0; pos < fragData.Size(); ++pos){
		scoresVec[pos] = prior * bgmodel[fragData[pos]];
	}
	return scoresVec;
}


