#include "functions_export_rcpp.hpp"



Rcpp::List calcStartCoverProbs_cpp(const Rcpp::IntegerVector& fragIDs,     // vector with unique fragment IDs
                                   const Rcpp::IntegerVector& fragPos,     // vector with positions within each fragment, 1 - based!
                                   const Rcpp::IntegerVector& protectVec,  // vector with protection data, 0 - accessible; 1 - protected
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


	// set parameters
	if(_VERBOSE_){
		Rcpp::Rcout<<"Creating PARAMS object..."<<endl;
	}
	parameters params(bgprotectprob,
                   bgprior);


	// create object with background/footprint models
	if(_VERBOSE_){
		Rcpp::Rcout<<"Creating footprint models object..."<<endl;
	}
	DNAbind_obj_vector ftp_models(binding_models,
                               params);



	// create object with SMF data
	if(_VERBOSE_){
		Rcpp::Rcout<<"Creating SEQUENCES object..."<<endl;
	}
	SMFdataset SMFdata(fragIDs,
                    fragPos,
                    protectVec,
                    ftp_models.maxwmlen);

	if(_VERBOSE_){
		Rcpp::Rcout<<"Creating Predict object..."<<endl;
	}

	// 	Predict predict(SMFdata,
	//                  ftp_models,
	//                  params);


	// run prediction
	if(_VERBOSE_){
		Rcpp::Rcout<<"Calculating posterior probabilities..."<<endl;
	}
	Predict predict;
	Rcpp::List outList = predict.calcStartCoverProbs(SMFdata,
                                                  ftp_models,
                                                  params,
                                                  report_prediction_in_flanks_,
                                                  Ncpu_);
	// clear
	SMFdata.clear();
	ftp_models.clear();

	return outList;
}


Rcpp::NumericMatrix count_spacing_freq_cpp(const Rcpp::IntegerVector& fragIDs,     // vector with unique fragment IDs
                                  const Rcpp::IntegerVector& fragPos,     // vector with positions within each fragment, 0 - based!
                                  const Rcpp::IntegerVector& protectVec,  // vector with protection data, 0 - accessible; 1 - protected
                                  const Rcpp::IntegerVector& maxspacing,
                                  const Rcpp::NumericVector& Ncpu,
                                  const Rcpp::LogicalVector& verbose){

	//set verbose
	extern bool _VERBOSE_;
	_VERBOSE_ = Rcpp::as<bool >(verbose);

	int Ncpu_ = Rcpp::as<int >(Ncpu);
#ifndef _OPENMP
	Rcpp::Rcout<<"nomeR was compiled without OpenMP. ncpu does not have effect.\n";
#endif


	int maxspacing_ = Rcpp::as<int >(maxspacing);

	SMFdataset SMFdata(fragIDs,
                    fragPos,
                    protectVec,
                    0);

	vector<vector<uint64_t> > freq_mat = SMFdata.count_freq_for_spacings(maxspacing_,
                                                                      Ncpu_);
	Rcpp::NumericMatrix ctable_out(maxspacing_,4);
	// Set row and column names
	Rcpp::colnames(ctable_out) = Rcpp::CharacterVector::create("N00", "N01", "N10", "N11");
	Rcpp::CharacterVector rnames(maxspacing_);
	for (int i = 0; i < maxspacing_; ++i)
		rnames[i] = to_string(i + 1);
	Rcpp::rownames(ctable_out) = rnames;
	for(size_t s=0; s<freq_mat.size(); ++s){
		ctable_out(s,0) = freq_mat[s][0];
		ctable_out(s,1) = freq_mat[s][1];
		ctable_out(s,2) = freq_mat[s][2];
		ctable_out(s,3) = freq_mat[s][3];
	}


	// //vector<uint64_t > spacings;
	// vector<uint64_t > freq00;
	// vector<uint64_t > freq01;
	// vector<uint64_t > freq10;
	// vector<uint64_t > freq11;
	//
	// for(int s=0; s<freq_mat.size(); ++s){
	// 	spacings.push_back(s + 1);
	// 	freq00.push_back(freq_mat[s][0]);
	// 	freq01.push_back(freq_mat[s][1]);
	// 	freq10.push_back(freq_mat[s][2]);
	// 	freq11.push_back(freq_mat[s][3]);
	// }
	//
	// Rcpp::List export_list;
	//
	// export_list.push_back(Rcpp::wrap(spacings),"S");
	// export_list.push_back(Rcpp::wrap(freq00),"N00");
	// export_list.push_back(Rcpp::wrap(freq01),"N01");
	// export_list.push_back(Rcpp::wrap(freq10),"N10");
	// export_list.push_back(Rcpp::wrap(freq11),"N11");

	SMFdata.clear();
	return ctable_out;

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


	for(size_t i = 0; i < p_joint.size(); ++i){
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




