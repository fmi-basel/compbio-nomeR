#include "parameters-class.hpp"

parameters::parameters(){

}


parameters::parameters(double _bgcoverprob,
           double _bgprior){
  setParams(_bgcoverprob,
            _bgprior);
}

parameters::~parameters()
{

}


void parameters::setParams(double _bgcoverprob,
                           double _bgprior){
  printoutonly = "All";
  bgcoverprob = _bgcoverprob;
  bgprior = _bgprior;

}

parameters & parameters::operator = (const parameters & other){
	if (this != &other){
    printoutonly = other.printoutonly;

    bgcoverprob = other.bgcoverprob;
    bgprior = other.bgprior;
	}
	return *this;
}

void parameters::print()
{
  Rcpp::Rcout <<"printoutonly\t"<<printoutonly<<endl;
  Rcpp::Rcout <<"bgcoverprob\t"<<bgcoverprob<<endl;
  Rcpp::Rcout <<"bgprior\t"<<bgprior<<endl;
}
void parameters::clear(){
  printoutonly = "All";// will be always "All" in R wrapper
  bgcoverprob = 0;
  bgprior= 0;

}
