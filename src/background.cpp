#include "background.h"

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
  cout << "Background parameters:\n";
  cout << "Background cover probability = "<<bgcoverprob<<endl;
  cout << "Background prior = "<<prior<<endl;
  
}
void Background::print_normalized() const{
  print();

}


double Background::get_score(int seq, int position) const{
   
  extern NOMeSeqData SEQUENCES;
  if(seq<0 || seq>=SEQUENCES.Size()){
	cerr<<"Wm::get_score: Index of sequence is out of range: "<<seq<<endl;
	exit(1);
  }

  double score = 1;
  for(int i = position;i < position + len;++i){
    if(i >= 0 && i < SEQUENCES[seq].Size())
      score *= bgmodel[SEQUENCES[seq][i]];
  }

  return prior * score;
  
}
