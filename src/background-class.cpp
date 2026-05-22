#include "background-class.hpp"

Background::Background(const parameters& params){
  classname = "background";
  name = "background";
  group = "background";
  len = 1;
  prior = params.bgprior;
  initialprior = prior;

  bgcoverprob = params.bgcoverprob;

  bgmodel.push_back(1-bgcoverprob); // [0]: emission for fully accessible position (mod_prob = 1)
  bgmodel.push_back(bgcoverprob);   // [1]: emission for fully protected position  (mod_prob = 0)
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
  for(int i = position; i < position + len; ++i){
    if(i >= 0 && i < SEQUENCES[seq].Size()){
      double p = SEQUENCES[seq][i];
      score *= (p < 0.0) ? 1.0 : (p * bgmodel[0] + (1.0 - p) * bgmodel[1]);
    }
  }
  return prior * score;
}

vector<double > Background::get_seq_scores_vec(const fragProtectData& fragData) const{
	vector<double > scoresVec(fragData.Size(), prior);
	for(int pos = 0; pos < (int)fragData.Size(); ++pos){
		double p = fragData[pos];
		scoresVec[pos] = prior * ((p < 0.0) ? 1.0 : (p * bgmodel[0] + (1.0 - p) * bgmodel[1]));
	}
	return scoresVec;
}

void Background::get_seq_scores_vec(const fragProtectData& fragData, vector<double>& out) const{
	out.assign(fragData.Size(), prior);
	for(int pos = 0; pos < (int)fragData.Size(); ++pos){
		double p = fragData[pos];
		out[pos] = prior * ((p < 0.0) ? 1.0 : (p * bgmodel[0] + (1.0 - p) * bgmodel[1]));
	}
}


