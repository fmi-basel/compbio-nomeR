#include "binding_object_model-class.hpp"

binding_object_model::binding_object_model(const vector<double > &_protect_prob,
                     const double _prior,
                     const string _name,
                     const string _group){

  classname = "wm";
  prior = _prior;
  initialprior = _prior;
  name = _name;
  group = _group;

  for(int pos = 0; pos < (int)_protect_prob.size(); ++pos){
    vector<double > tmp;
    tmp.push_back(1 - _protect_prob[pos]); // emit for fully accessible position (mod_prob = 1)
    tmp.push_back(_protect_prob[pos]);     // emit for fully protected position  (mod_prob = 0)
    mat.push_back(tmp);
  }
  len = mat.size();
  normalize();

}

binding_object_model::~binding_object_model() {

}

void binding_object_model::print() const{
  Rcpp::Rcout <<"//\n";
  Rcpp::Rcout <<"NA\t"<< name <<"\n";
  Rcpp::Rcout <<"GROUP\t"<< group <<"\n";

  Rcpp::Rcout <<"PRIOR\t"<< prior <<"\n";
  Rcpp::Rcout <<"POS\tUNPROT\tPROT\n";
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
  Rcpp::Rcout <<"GROUP\t"<< group <<"\n";
  //Rcpp::Rcout <<"Orientation\t"<<orientation<<endl;
  Rcpp::Rcout <<"PRIOR\t"<< prior <<"\n";
  Rcpp::Rcout <<"POS\tUNPROT\tPROT\n";
  for(int i=0;i<normmat.size();i++){
    Rcpp::Rcout <<i+1;
    for(int j=0;j<normmat[i].size();j++)
      Rcpp::Rcout << "\t"<< normmat[i][j];
    Rcpp::Rcout <<"\n";
  }
  Rcpp::Rcout <<"//\n";
}


double binding_object_model::get_score(const SMFdataset& SEQUENCES,
                                       int seq,
                                       int position) const{


  if(seq<0 || seq>=SEQUENCES.Size()){
    Rcpp::stop("binding_object_model::get_score: Index of sequence is out of range: ");
	}

  double score = 1;
  for(int i = position; i < position + len; ++i){
    if(i >= 0 && i < SEQUENCES[seq].Size()){
      double p = SEQUENCES[seq][i];
      int mpos = i - position;
      score *= (p < 0.0) ? 1.0 : (p * normmat[mpos][0] + (1.0 - p) * normmat[mpos][1]);
    }
  }
  return prior * score;
}

vector<double > binding_object_model::get_seq_scores_vec(const fragProtectData& fragData) const{
	vector<double > scoresVec(fragData.Size(), prior);
	// score at position 0: full product over all model positions
	double score = 1.0;
	for(int i = 0; i < len; ++i){
		if(i < (int)fragData.Size()){
			double p = fragData[i];
			score *= (p < 0.0) ? 1.0 : (p * normmat[i][0] + (1.0 - p) * normmat[i][1]);
		}
	}
	scoresVec[0] = prior * score;
	// sliding window O(n): score(pos) = score(pos-1) * emit(new_last) / emit(old_first)
	// exact for uniform models (normmat[0] == normmat[len-1]), which is standard in SMF
	for(int pos = 1; pos <= (int)fragData.Size() - len; ++pos){
		double p_new = fragData[pos + len - 1];
		double e_new = (p_new < 0.0) ? 1.0 : (p_new * normmat[len-1][0] + (1.0 - p_new) * normmat[len-1][1]);
		double p_old = fragData[pos - 1];
		double e_old = (p_old < 0.0) ? 1.0 : (p_old * normmat[0][0] + (1.0 - p_old) * normmat[0][1]);
		scoresVec[pos] = scoresVec[pos - 1] * e_new / e_old;
	}
	return scoresVec;
}

void binding_object_model::get_seq_scores_vec(const fragProtectData& fragData, vector<double>& out) const{
	out.assign(fragData.Size(), prior);
	double score = 1.0;
	for(int i = 0; i < len; ++i){
		if(i < (int)fragData.Size()){
			double p = fragData[i];
			score *= (p < 0.0) ? 1.0 : (p * normmat[i][0] + (1.0 - p) * normmat[i][1]);
		}
	}
	out[0] = prior * score;
	for(int pos = 1; pos <= (int)fragData.Size() - len; ++pos){
		double p_new = fragData[pos + len - 1];
		double e_new = (p_new < 0.0) ? 1.0 : (p_new * normmat[len-1][0] + (1.0 - p_new) * normmat[len-1][1]);
		double p_old = fragData[pos - 1];
		double e_old = (p_old < 0.0) ? 1.0 : (p_old * normmat[0][0] + (1.0 - p_old) * normmat[0][1]);
		out[pos] = out[pos - 1] * e_new / e_old;
	}
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
