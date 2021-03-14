#include "DNAbindobj_vector.h"

DNAbind_obj_vector::DNAbind_obj_vector(){
  
}

DNAbind_obj_vector::~DNAbind_obj_vector(){
  
  for(int wm=0;wm<objvector.size();wm++){
  	if(objvector[wm] != NULL){
  		delete objvector[wm];
  	}
  }
}

DNAbinding_object * DNAbind_obj_vector::operator [](int i){
	if(i>objvector.size()-1 || i<0){
	  Rcpp::stop("DNAbind_obj_vector::operator[](int i):  The index is out of range");
		// Rcpp::Rcerr<<"DNAbind_obj_vector::operator[](int i):  The index is out of range: "<<i<<"\n";
		// exit(1);
	}
	return objvector[i];
}

int DNAbind_obj_vector::Size(){
	return size;
}


int DNAbind_obj_vector::create(Rcpp::List _bind_objs,
                               parameters &params){
  numberofobjects=0;
  maxwmlen = 1;
  //cout<<"#### In DNAbind_obj_vector create ####"<<endl;
  //cout<<"List size is "<<_bind_objs.size()<<endl;
  for(int wm=0; wm < _bind_objs.size(); ++wm){
    Rcpp::List pinf = Rcpp::as<Rcpp::List >(_bind_objs[wm]);
    //cout<<"Index="<<wm<<" Data list size "<<pinf.size()<<endl;
    if(!pinf.containsElementNamed("PROTECT_PROB")){
      Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element PROTECT_PROB\n";
      return(0);
    }
    if(!pinf.containsElementNamed("PRIOR")){
      Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element PRIOR\n";
      return(0);
    }
    if(!pinf.containsElementNamed("NAME")){
      Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element NAME\n";
      return(0);
    }
    
    vector<double > prot_prob = Rcpp::as<vector<double > >(pinf["PROTECT_PROB"]);
    double prior = Rcpp::as<double >(pinf["PRIOR"]);
    string name = Rcpp::as<string >(pinf["NAME"]);

    DNAbinding_object *newobj=new binding_object_model(prot_prob,
                                                       prior,
                                                       name);
    objvector.push_back(newobj);
    names.push_back(newobj->name);
    priors.push_back(newobj->prior);
    
    if(newobj->len>maxwmlen)
      maxwmlen = newobj->len;
    
    vector<int > tmp(1,0);

    tmp[0] = objvector.size() - 1;
    names2index.push_back(tmp);
    
  }
  
  
  //initialise backgound model
  
  Background *bg = new Background(params);
  objvector.push_back(bg);
  names.push_back(bg->name);
  priors.push_back(bg->prior);
  
  vector<int > tmp(1,objvector.size() - 1);
  names2index.push_back(tmp);
  
  numberofobjects = names.size();
  size = objvector.size();
  
  
  // normalize priors so that they sum up to 1
  
  double priorsum=0;
  for(int i=0;i<objvector.size();++i){
    priorsum += objvector[i]->prior;// * objvector[i]->len;
  }
  
  for(int i=0;i<objvector.size();++i){
    objvector[i]->prior = (objvector[i]->prior)/priorsum;
    priors[i] = objvector[i]->prior;
    //objvector[i]->print_normalized();
  }
  
  return objvector.size();
  
  
}

void DNAbind_obj_vector::clear(){
  numberofobjects=0;
  maxwmlen=0;
  names.clear();
  names2index.clear();
  priors.clear();
  
  for(int wm=0;wm<objvector.size();wm++){
    if(objvector[wm] != NULL){
      delete objvector[wm];
    }
  }
  objvector.clear();
  size=0;
  
  
  
}


void DNAbind_obj_vector::print(){
  
  for(int wm = 0;wm < objvector.size(); ++wm){
    objvector[wm]->print_normalized();
  }
  
}





// Function to calculate theoretical joint probabiltiies to observe 00,01,10,11 at distance S

vector<vector<double > > DNAbind_obj_vector::calc_theor_joint_prob(vector<double > ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector
                                               // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
                                               double bg_protect_prob,
                                               double footprint_protect_prob,
                                               int max_spacing){
  
  // define vector with ftp lengths
  vector<int > ftp_lengths;
  double R_const=0;
  for(int ftp=0;ftp < ftp_cover_priors.size();++ftp){
    ftp_lengths.push_back(ftp + 1);
    R_const += ftp_cover_priors[ftp]/(ftp + 1);
  }
  
  if(R_const <= 0){
    Rcpp::stop("DNAbind_obj_vector::calc_theor_joint_prob: ERROR! Value of R_const is negative or zero.");
  }
  // calculate start priors for ftps
  vector<double > ftp_start_priors;
  for(int ftp=0; ftp < ftp_cover_priors.size();++ftp){
    ftp_start_priors.push_back((ftp_cover_priors[ftp]/ftp_lengths[ftp])/R_const);
  }
  
  // calculate forward partition sum
  vector<double > FW(max_spacing + 1,0); // here FW[0] is intial condition for FW
  vector<double > ProbBg(max_spacing + 1,0);
  vector<double > sigma_cumul(max_spacing + 1,0);
  FW[0] = 1;
  ProbBg[0] = 0;
  sigma_cumul[0] = 0;
  for(int d=1; d <= max_spacing; ++d){
    FW[d] = 0;
    for(int ftp=0; ftp < ftp_start_priors.size(); ++ftp){
      if(d - ftp_lengths[ftp] >=0){
        FW[d] += ftp_start_priors[ftp] * FW[d - ftp_lengths[ftp]];
      }
    }
    
    // set bg prob
    ProbBg[d] = ftp_start_priors[0] * FW[d-1];
    sigma_cumul[d] = sigma_cumul[d-1] + ProbBg[d];
  }
  
  // set emission probs vectors for convenience
  double beta1 = bg_protect_prob;
  double beta0 = 1-bg_protect_prob;
  double alpha1 = footprint_protect_prob;
  double alpha0 = 1-footprint_protect_prob;
  vector<double > alpha1_vec;
  vector<double > alpha0_vec;
  alpha1_vec.push_back(beta1);
  alpha0_vec.push_back(beta0);
  for(int ftp=1; ftp < ftp_start_priors.size(); ++ftp){
    alpha1_vec.push_back(alpha1);
    alpha0_vec.push_back(alpha0);
  }
  
  // calculate joint theoretical probabilities
  vector<vector<double > > Pjoint;
  // fill at S=1, i.e probs of 1 and 0. no distance
  vector<double > tmp(5,0);
  tmp[0] = 1;
  tmp[1] = alpha0 + (beta0 - alpha0) * ftp_cover_priors[0];
  tmp[2] = 0;
  tmp[3] = 0;
  tmp[4] = alpha1 + (beta1 - alpha1) * ftp_cover_priors[0];
  Pjoint.push_back(tmp);
  for(int S = 2; S <= max_spacing; ++S){
    vector<double > prob_joint(5,0); // This correspond to Spacing, P(0,0|S), P(0,1|S), P(1,0|S), P(1,1|S)
    double sigma_sum0 = 0;
    double sigma_sum1 = 0;
    for(int ftp=0; ftp < ftp_start_priors.size(); ++ftp){
      int ftp_len = ftp_lengths[ftp];
      if(S - ftp_len - 1 >0){
        sigma_sum0 += ftp_start_priors[ftp] * alpha0_vec[ftp] * sigma_cumul[S - ftp_len - 1];
        sigma_sum1 += ftp_start_priors[ftp] * alpha1_vec[ftp] * sigma_cumul[S - ftp_len - 1];
      }
    }
    
    prob_joint[0] = S;
    // calc probs
    
    // P(0,0|S)
    prob_joint[1] = alpha0 * alpha0 + R_const * ftp_start_priors[0] * alpha0 * (beta0 - alpha0) +
      R_const * (beta0 - alpha0) * ( (alpha0 + (beta0 - alpha0) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum0);
    // P(0,1|S)
    prob_joint[2] = alpha1 * alpha0 + R_const * ftp_start_priors[0] * alpha1 * (beta0 - alpha0) + 
      R_const * (beta1 - alpha1) * ( (alpha0 + (beta0 - alpha0) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum0);
    // P(1,0|S)
    prob_joint[3] = alpha0 * alpha1 + R_const * ftp_start_priors[0] * alpha0 * (beta1 - alpha1) + 
      R_const * (beta0 - alpha0) * ( (alpha1 + (beta1 - alpha1) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum1);
    // P(1,1|S)
    prob_joint[4] = alpha1 * alpha1 + R_const * ftp_start_priors[0] * alpha1 * (beta1 - alpha1) + 
      R_const * (beta1 - alpha1) * ( (alpha1 + (beta1 - alpha1) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum1);
    
    
    Pjoint.push_back(prob_joint);
  }

  return(Pjoint);
  
}

