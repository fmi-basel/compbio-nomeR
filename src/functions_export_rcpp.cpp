#include "functions_export_rcpp.hpp"



Rcpp::List run_cpp_nomeR(const Rcpp::List& data,
                         const Rcpp::CharacterVector& fragnames,
                         const Rcpp::List& binding_models,
                         const Rcpp::NumericVector& bgprotectprob,
                         const Rcpp::NumericVector& bgprior,
                         const Rcpp::LogicalVector& report_prediction_in_flanks,
                         const Rcpp::NumericVector& Ncpu,
                         const Rcpp::LogicalVector& verbose
) {

	//set verbose
	extern bool _VERBOSE_;
	_VERBOSE_ = Rcpp::as<bool >(verbose);

	// set report_prediction_in_flanks
	bool report_prediction_in_flanks_ = Rcpp::as<bool >(report_prediction_in_flanks);

	int Ncpu_ = Rcpp::as<int >(Ncpu);
#ifndef _OPENMP
	Rcpp::Rcout<<"nomeR was compiled without OpenMP. ncpu does not have effect.\n";
#endif



	if(_VERBOSE_){
		Rcpp::Rcout<<"Creating Predict object..."<<endl;
	}
	Predict predict(data,
                 fragnames,
                 binding_models,
                 bgprotectprob,
                 bgprior);

	// run prediction
	if(_VERBOSE_){
		Rcpp::Rcout<<"Calculating posterior binding probabilities..."<<endl;
	}

	Rcpp::List output_data;
	if(predict.Run(Ncpu_)){
		//predict.Run(Ncpu_);

		if(_VERBOSE_){
			Rcpp::Rcout<<"Running predict.getStartProbDF()..."<<endl;
		}
		Rcpp::List startProbdf = predict.getStartProbDF(report_prediction_in_flanks_);
		if(_VERBOSE_){
			Rcpp::Rcout<<"Running predict.getCoverProbDF()..."<<endl;
		}
		Rcpp::List coverProb = predict.getCoverProbDF();

		output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = startProbdf,
                                    Rcpp::Named("COVER_PROB") = coverProb);


// 		if(_VERBOSE_){
// 			Rcpp::Rcout<<"Running predict.getGenomeSummaryDF()..."<<endl;
// 		}
// 		Rcpp::List genSummary = predict.getGenomeSummaryDF();
// 		output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = startProbdf,
//                                     Rcpp::Named("COVER_PROB") = coverProb,
//                                     Rcpp::Named("SUMMARY") = genSummary);
	} else {
		output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = R_NilValue,
                                    Rcpp::Named("COVER_PROB") = R_NilValue);

	}

	_VERBOSE_ = 0;
	return output_data;

}


Rcpp::List count_spacing_freq_cpp(const Rcpp::List& data,
                                  const Rcpp::CharacterVector& fragnames,
                                  const Rcpp::IntegerVector& maxspacing,
                                  const Rcpp::IntegerVector& maxwmlen){

	int maxwmlen_ = Rcpp::as<int >(maxwmlen);
	int maxspacing_ = Rcpp::as<int >(maxspacing);

	NOMeSeqData nome_data(data,
                       fragnames,
                       maxwmlen_);
	//nome_data.create(data,maxwmlen_);
	vector<vector<int> > freq_mat = nome_data.count_freq_for_spacings(maxspacing_);
	vector<int > spacings;
	vector<int > freq00;
	vector<int > freq01;
	vector<int > freq02;
	vector<int > freq10;
	vector<int > freq11;
	vector<int > freq12;

	vector<int > freq20;
	vector<int > freq21;
	vector<int > freq22;
	for(int i=0; i<freq_mat.size(); ++i){
		spacings.push_back(freq_mat[i][0]);
		freq00.push_back(freq_mat[i][1]);
		freq01.push_back(freq_mat[i][2]);
		freq02.push_back(freq_mat[i][3]);
		freq10.push_back(freq_mat[i][4]);
		freq11.push_back(freq_mat[i][5]);
		freq12.push_back(freq_mat[i][6]);
		freq20.push_back(freq_mat[i][7]);
		freq21.push_back(freq_mat[i][8]);
		freq22.push_back(freq_mat[i][9]);

	}

	Rcpp::List export_list;

	export_list.push_back(Rcpp::wrap(spacings),"S");
	export_list.push_back(Rcpp::wrap(freq00),"N00");
	export_list.push_back(Rcpp::wrap(freq01),"N01");
	export_list.push_back(Rcpp::wrap(freq10),"N10");
	export_list.push_back(Rcpp::wrap(freq11),"N11");
	export_list.push_back(Rcpp::wrap(freq02),"N0na");
	export_list.push_back(Rcpp::wrap(freq12),"N1na");
	export_list.push_back(Rcpp::wrap(freq20),"Nna0");
	export_list.push_back(Rcpp::wrap(freq21),"Nna1");
	export_list.push_back(Rcpp::wrap(freq22),"Nnana");


	nome_data.clear();
	return export_list;

}


Rcpp::List calculate_theor_joint_prob_cpp(const Rcpp::NumericVector& ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                          const Rcpp::NumericVector& bg_protect_prob,
                                          const Rcpp::NumericVector& footprint_protect_prob,
                                          const Rcpp::IntegerVector& max_spacing){

	vector<double > ftp_cover_priors_ = Rcpp::as<vector<double > >(ftp_cover_priors);
	double bg_protect_prob_ = Rcpp::as<double >(bg_protect_prob);
	double footprint_protect_prob_ = Rcpp::as<double >(footprint_protect_prob);
	int max_spacing_ = Rcpp::as<int >(max_spacing);


	DNAbind_obj_vector ftp_array;
	vector<vector<double > > p_joint = ftp_array.calc_theor_joint_prob(ftp_cover_priors_,
                                                                    bg_protect_prob_,
                                                                    footprint_protect_prob_,
                                                                    max_spacing_);

	// create Rcpp object

	vector<int > spacings;
	vector<double > p00;
	vector<double > p01;
	vector<double > p10;
	vector<double > p11;


	for(int i = 0; i < p_joint.size(); ++i){
		//Rcpp::Rcout<<"S="<<p_joint[i][0]<<"; "<<p_joint[i][1]<<"; "<<p_joint[i][2]<<"; "<<p_joint[i][3]<<"; "<<p_joint[i][4]<<endl;
		spacings.push_back(p_joint[i][0]);
		p00.push_back(p_joint[i][1]);
		p01.push_back(p_joint[i][2]);
		p10.push_back(p_joint[i][3]);
		p11.push_back(p_joint[i][4]);
	}

	Rcpp::List export_list;
	export_list.push_back(Rcpp::wrap(spacings),"S");
	export_list.push_back(Rcpp::wrap(p00),"P00");
	export_list.push_back(Rcpp::wrap(p01),"P01");
	export_list.push_back(Rcpp::wrap(p10),"P10");
	export_list.push_back(Rcpp::wrap(p11),"P11");


	ftp_array.clear();
	return export_list;
}




