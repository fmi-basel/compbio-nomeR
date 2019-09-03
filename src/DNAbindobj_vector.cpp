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


