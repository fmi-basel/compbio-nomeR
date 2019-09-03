#include "parameters.h"

parameters::parameters(){

}


parameters::parameters(double _bgcoverprob,
           double _bgprior,
           double _bound_fit_tol,
           double _bound_min_fderiv_val,
           int _bound_max_steps,
           double _priorEM_fit_tol,
           int _priorEM_max_steps){
  setParams(_bgcoverprob,
            _bgprior,
            _bound_fit_tol,
            _bound_min_fderiv_val,
            _bound_max_steps,
            _priorEM_fit_tol,
            _priorEM_max_steps);
}
  
parameters::~parameters()
{

}


void parameters::setParams(double _bgcoverprob,
                           double _bgprior,
                           double _bound_fit_tol,
                           double _bound_min_fderiv_val,
                           int _bound_max_steps,
                           double _priorEM_fit_tol,
                           int _priorEM_max_steps){
  printoutonly = "All";
  BOUND_FIT_TOLERANCE = _bound_fit_tol;
  BOUND_MIN_FDERIV_VAL = _bound_min_fderiv_val;
  BOUND_MAX_STEPS = _bound_max_steps;
  
  PRIOREM_FIT_TOLERANCE = _priorEM_fit_tol;
  PRIOREM_MAX_STEPS = _priorEM_max_steps;
  
  bgcoverprob = _bgcoverprob;
  bgprior = _bgprior;
  
}

parameters & parameters::operator = (const parameters & other){
	if (this != &other){
    printoutonly = other.printoutonly;

    bgcoverprob = other.bgcoverprob;
    bgprior = other.bgprior;

    BOUND_FIT_TOLERANCE = other.BOUND_FIT_TOLERANCE;
    BOUND_MAX_STEPS = other.BOUND_MAX_STEPS;
    BOUND_MIN_FDERIV_VAL = other.BOUND_MIN_FDERIV_VAL;

    PRIOREM_FIT_TOLERANCE = other.PRIOREM_FIT_TOLERANCE;
    PRIOREM_MAX_STEPS = other.PRIOREM_MAX_STEPS;

	}
	return *this;
}

void parameters::print()
{
  Rcpp::Rcout <<"printoutonly\t"<<printoutonly<<endl;
  Rcpp::Rcout <<"bgcoverprob\t"<<bgcoverprob<<endl;
  Rcpp::Rcout <<"bgprior\t"<<bgprior<<endl;
  Rcpp::Rcout <<"BOUND_FIT_TOLERANCE\t"<<BOUND_FIT_TOLERANCE<<endl;
  Rcpp::Rcout <<"BOUND_MAX_STEPS\t"<<BOUND_MAX_STEPS<<endl;
  Rcpp::Rcout <<"BOUND_MIN_FDERIV_VAL\t"<<BOUND_MIN_FDERIV_VAL<<endl;
  Rcpp::Rcout <<"PRIOREM_FIT_TOLERANCE\t"<<PRIOREM_FIT_TOLERANCE<<endl;
  Rcpp::Rcout <<"PRIOREM_MAX_STEPS\t"<<PRIOREM_MAX_STEPS<<endl;


}
void parameters::clear(){
  printoutonly = "All";// will be always "All" in R wrapper
  BOUND_FIT_TOLERANCE = 0;
  BOUND_MIN_FDERIV_VAL = 0;
  BOUND_MAX_STEPS = 0;
  
  PRIOREM_FIT_TOLERANCE = 0;
  PRIOREM_MAX_STEPS = 0;
  
  bgcoverprob = 0;
  bgprior= 0;
  
}