#include "predict-class.hpp"

//Predict::Predict(){}
Predict::~Predict(){
	
}

Predict::Predict(const SMFdataset& refSmfData,
                 const DNAbind_obj_vector& refFtp_models,
                 const parameters& refParams)
	: SEQUENCES(refSmfData),
   BINDING_OBJECTS(refFtp_models),
   PARAMS(refParams)
{
	Create();
}


bool Predict::Create()
{
	// verbosity
	extern bool _VERBOSE_;
	
	// Initialize print_indexes
	if(PARAMS.printoutonly != "All"){
		Tokenize(PARAMS.printoutonly,print_names," \t,;{}");
		
		for(int i=0;i<print_names.size();++i){
			
			for(int wm=0;wm<BINDING_OBJECTS.Size();++wm){
				if(BINDING_OBJECTS[wm] -> name == print_names[i]){
					print_indexes.push_back(wm);
				}
			}
		}
		if(print_indexes.size() == 0){
			Rcpp::stop("Forward_Backward_algorithm:: ERROR! Incorrect value of parameter printoutonly. No objects to calculate were found.\n");
		}
	}
	else{
		for(int wm = 0;wm < BINDING_OBJECTS.Size(); ++wm){
			print_indexes.push_back(wm);
			vector<string >::iterator it;
			it = find(print_names.begin(),print_names.end(),BINDING_OBJECTS[wm]->name);
			if(it == print_names.end())
				print_names.push_back(BINDING_OBJECTS[wm]->name);
		}
	}
	
	// define maps indexes2names and names2indexes
	// names2indexes - this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
	for(int i=0;i<print_names.size();++i){
		vector<int > tmp;
		for(int wm=0;wm<print_indexes.size();++wm){
			if(BINDING_OBJECTS[print_indexes[wm]]->name == print_names[i])
				tmp.push_back(print_indexes[wm]);
		}
		names2indexes.push_back(tmp);
	}
	
	// define names2indicesinprobarray
	for(int i=0;i<print_names.size();++i){
		vector<int > tmp;
		for(int wm=0;wm<print_indexes.size();++wm){
			if(BINDING_OBJECTS[print_indexes[wm]]->name == print_names[i])
				tmp.push_back(wm);
		}
		names2indicesinprobarray.push_back(tmp);
		
	}
	return 1;
}



// method that calculates start and cover probabilities and returns a Rcpp::List with calculated data.

Rcpp::List Predict::calcStartCoverProbs(bool report_prediction_in_flanks,
                                        int ncpu){
	extern bool _VERBOSE_;
	
	int numberofobjects = BINDING_OBJECTS.Size(); // number of footprint models including background
	int maxwmlen = BINDING_OBJECTS.maxwmlen; // maximum size of footprint model;
	int seq = 0;
	// initial value for partition sums.
	// when footprint priors are normalized, i.e. sum of all priors is 1, then initial values for partition sums is always 1.
	// we normalize the priors, therefore we set value to 1.
	double part_init = 1;
	
	// allocate natice C++ vectors for output probabilities
	vector<int32_t > tmp_vec;
	vector<vector<int32_t >> startOutFragIDs(SEQUENCES.Size(),tmp_vec); // vector of vectors with fragment IDs. one per seq
	vector<vector<int32_t >> startOutFragPos(SEQUENCES.Size(),tmp_vec); // positions within fragments
	
	vector<vector<int32_t >> coverOutFragIDs(SEQUENCES.Size(),tmp_vec); // vector with fragment IDs as was passed from the R side
	vector<vector<int32_t >> coverOutFragPos(SEQUENCES.Size(),tmp_vec); // positions within fragments
	
	vector<vector<vector<double >>> startOutProbs; // vectors of size print_names, i.e. for each footprint. per each seq
	vector<vector<vector<double >>> coverOutProbs;
	
	for(seq = 0; seq < SEQUENCES.Size(); ++seq){
		vector<vector<double >> tmpst;
		vector<vector<double >> tmpcv;
		for(int i=0; i<print_names.size(); ++i){
			vector<double > startTmp;
			vector<double > coverTmp;
			tmpst.push_back(startTmp);
			tmpcv.push_back(coverTmp);
		}
		startOutProbs.push_back(tmpst);
		coverOutProbs.push_back(tmpcv);
	}
	// all vectors above, i.e. 
	// startOutFragIDs, startOutFragPos, startOutProbs[ftp]
	// and coverOutFragIDs, coverOutFragPos, coverOutProbs[ftp]
	// will be appended during calculation of probabilities for each molecule.
	// update: appending is not thread-safe!!! work with local vectors 
	
#ifdef _OPENMP
	omp_set_nested(true);
	omp_set_num_threads(ncpu);
	if(_VERBOSE_)
		Rcpp::Rcout<<"Running prediction with "<<omp_get_max_threads()<<" cpu."<<endl;
#endif
	
#pragma omp parallel private(seq)
{
	
#pragma omp for schedule(dynamic)
	for(seq = 0; seq < SEQUENCES.Size(); ++seq){
		int seqlength = SEQUENCES[seq].Size();
		
		// allocate memory for:
		// F -  forward parition sum
		// R - backward partition sum
		// Prob - probability of footprint ends at position pos
		
		// this is how it was allocated previously
		// Allocate memory for F, R and Prob
		// int maxwmlen = BINDING_OBJECTS.maxwmlen;
		// for(int seq=0; seq < SEQUENCES.Size();++seq){
		//   vector<double > tmp(SEQUENCES[seq].Size() + maxwmlen + 1,1);
		//   F.push_back(tmp);
		//   tmp.resize(SEQUENCES[seq].Size() + 2,1);
		//   R.push_back(tmp);
		//   vector<double > tmp2(SEQUENCES[seq].Size() + 1, 0);
		//   vector<vector<double > > tmp1(print_indexes.size(),tmp2);
		//   Prob.push_back(tmp1);
		// 
		// }
		vector<double > F(SEQUENCES[seq].Size() + maxwmlen + 1,1); // allocate memory for forward parition sum
		vector<double > R(SEQUENCES[seq].Size() + 2,1); // allocate memory for backward partition sum. it is shorter than F
		vector<double > probPerFtp(SEQUENCES[seq].Size() + 1, 0);
		vector<vector<double > > Prob(print_indexes.size(),probPerFtp);
		
		// calculate forward partition summ
		F[0] = part_init;
		vector<double > pf(numberofobjects,1);
		for(int pos=1;pos <= seqlength + BINDING_OBJECTS.maxwmlen; ++pos){
			double summ=0;
			for(int wm = 0; wm < numberofobjects; ++wm){
				pf[wm] = 1;
				int objlen = BINDING_OBJECTS[wm]->len;
				
				if(BINDING_OBJECTS[wm]->prior > 0){
					pf[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq, pos - objlen);
					for(int i = pos - objlen + 1; i <= pos - 1; ++i){
						
						if(i>=0){
							pf[wm] *= F[i];
						} else {
							pf[wm] *= part_init;
						}
					}
					
				}
				else{
					pf[wm] = 0;
				}
				summ += pf[wm];
			}
			
			F[pos] = 1/summ;
			if(pos<=seqlength){
				for(int wm=0;wm<print_indexes.size();++wm){
					Prob[wm][pos] = pf[print_indexes[wm]];
				}
			}
		}
		
		// calculate backward partition summ
		vector<double > pb(numberofobjects, 1);
		R[seqlength + 1] = part_init;
		for(int pos = seqlength; pos>=1; --pos){
			double summ = 0;
			for(int wm = 0; wm < numberofobjects; ++wm){
				pb[wm] = 1;
				int objlen = BINDING_OBJECTS[wm]->len;
				if(BINDING_OBJECTS[wm]->prior > 0){
					pb[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq,pos - 1);
					for(int i = pos+1; i <= pos+objlen-1; ++i){
						if(i<=seqlength + 1){
							pb[wm] *= R[i];
						} else{
							pb[wm] *= part_init;
						}
					}
					
				}
				else{
					pb[wm] = 0;
				}
				summ += pb[wm];
			}
			
			R[pos] = 1/summ;
		}
		// сalculate Z = Fn/Rn
		
		// calculate initial value for Z based on requirement that total coverage at L must be 1
		double zsumm =0;
		for(int wm = 0; wm < numberofobjects; ++wm){
			int objlen = BINDING_OBJECTS[wm]->len;
			double wmsumm = 0;
			for(int pos = seqlength; pos <= seqlength + objlen - 1; ++pos){
				double prod = 1;
				for(int j=pos - objlen + 1; j<=seqlength;++j){
					prod *= F[j];
				}
				wmsumm += pow(part_init,pos - seqlength) * prod;
			}
			zsumm += BINDING_OBJECTS[wm]->prior * wmsumm;
		}
		
		double z_init = 1/zsumm;
		
		R[seqlength + 1] = z_init;
		for(int pos = seqlength;pos >= 1;--pos){
			R[pos] = F[pos] * R[pos + 1]/R[pos];
		}
		// calculate start posteriors 
		for(int pos = 1;pos <= seqlength; ++pos){
			
			for(int wm=0;wm < print_indexes.size(); ++wm){
				int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
				if(pos + objlen - 1 <= seqlength)
					Prob[wm][pos] = Prob[wm][pos + objlen - 1] * F[pos + objlen -1] * R[pos+objlen];
				else
					Prob[wm][pos] = 0;
			}
		}
		
		
		// get output data structure for current sequence for start and cover probabilities
		// namely, append the vectors that we initilized before the for loop
		// i.e. startOutFragIDs, startOutFragPos, startOutProbs[ftp]
		// and coverOutFragIDs, coverOutFragPos, coverOutProbs[ftp]
		// for the current molecule
		
		
		
		// fill output vectors for START_PROB
		int firstDatPos = SEQUENCES[seq]._firstDatpos;
		int lastDatPos = SEQUENCES[seq]._lastDatpos;
		int spos = report_prediction_in_flanks ? 1 : firstDatPos;
		int lpos = lastDatPos;
		vector<int32_t > currSeqStartOutFragIDs;
		vector<int32_t > currSeqStartOutFragPos;
		vector<double > tmpstart;
		vector<vector<double >> currSeqStartOutProbs(print_names.size(),tmpstart);
		for(int position = spos; position <= lpos; ++position){
			//// 1. fill fragIDs and fragPos
			currSeqStartOutFragIDs.push_back(SEQUENCES[seq].Name());
			currSeqStartOutFragPos.push_back(position - firstDatPos + 1);
			//// 2. fill start probabilities for each footprint
			for(int i=0; i<print_names.size(); ++i){
				// summ across probablities associated with current ftp
				double totalprob=0;
				for(int j=0; j<names2indicesinprobarray[i].size(); ++j){
					totalprob += Prob[names2indicesinprobarray[i][j]][position];
				}
				currSeqStartOutProbs[i].push_back(totalprob);
			}
		}
		startOutFragIDs[seq] = move(currSeqStartOutFragIDs);
		startOutFragPos[seq] = move(currSeqStartOutFragPos);
		startOutProbs[seq] = move(currSeqStartOutProbs);
		
		
		// fill output vectors for COVER_PROB. Perhaps, this can be optimized by adding and subtracting start prob at end and beginning of footprint
		vector<int32_t > currSeqCoverOutFragIDs;
		vector<int32_t > currSeqCoverOutFragPos;
		vector<double > tmpcov;
		vector<vector<double >> currSeqCoverOutProbs(print_names.size(),tmpcov);
		
		for(int position = firstDatPos; position <= lastDatPos; ++position){
			//// 1. fill fragIDs and fragPos
			currSeqCoverOutFragIDs.push_back(SEQUENCES[seq].Name());
			currSeqCoverOutFragPos.push_back(position - firstDatPos + 1);
			//// 2. fill cover probabilities for each footprint
			for(int i=0; i<print_names.size(); ++i){
				double coverprob = 0;
				int objlen = BINDING_OBJECTS[names2indexes[i][0]]->len;
				for(int j=0; j<names2indicesinprobarray[i].size(); ++j){
					for(int p = max(position - objlen + 1,1); p <= position; ++p){
						coverprob += Prob[names2indicesinprobarray[i][j]][p];
					}
				}
				currSeqCoverOutProbs[i].push_back(coverprob);
			}
		}
		
		coverOutFragIDs[seq] = move(currSeqCoverOutFragIDs);
		coverOutFragPos[seq] = move(currSeqCoverOutFragPos);
		coverOutProbs[seq] = move(currSeqCoverOutProbs);
		
	}
	
} // end of omp parallel


// flatten nested vectors and create Rcpp::vectors for START_PROB
// memory for start probabilities
Rcpp::List RcppListStartOut; // this is a list of vectors
// 1st element: Rcpp::IntegerVector with fragment IDs as was passed from the R side
// 2nd element: Rcpp::IntegerVector with positions within fragments
// 3rd, 4th and so on: Rcpp::NumericVector with starting probabilities for ftp1, ftp2 and so on
// 1. Compute total length
size_t start_total_size = 0;
for (const auto& v : startOutFragIDs)
	start_total_size += v.size();
Rcpp::IntegerVector RcppStartOutFragIDs(start_total_size);
size_t offset = 0;
for (const auto& v : startOutFragIDs) {
	std::copy(v.begin(), v.end(), RcppStartOutFragIDs.begin() + offset);
	offset += v.size();
}
RcppListStartOut.push_back(RcppStartOutFragIDs,"seq");

Rcpp::IntegerVector RcppStartOutFragPos(start_total_size);
offset = 0;
for (const auto& v : startOutFragPos) {
	std::copy(v.begin(), v.end(), RcppStartOutFragPos.begin() + offset);
	offset += v.size();
}
RcppListStartOut.push_back(RcppStartOutFragPos,"pos");

// add flattened start probabilities for each footprint
for(int i=0; i<print_names.size(); ++i){
	Rcpp::NumericVector ftpStartProbs(start_total_size,NA_REAL);
	RcppListStartOut.push_back(ftpStartProbs,print_names[i]);
}
offset = 0;
for(int seq = 0; seq < startOutProbs.size(); ++seq){
	for(int i = 0; i < print_names.size(); ++i){
		Rcpp::NumericVector ftpProbVec = RcppListStartOut[i + 2]; // 0 - fragID, 1 - fragPos, 2 - ftp1, 3 - ftp2 etc.
		std::copy(startOutProbs[seq][i].begin(), startOutProbs[seq][i].end(), ftpProbVec.begin() + offset);
	}
	offset += startOutProbs[seq][0].size();
}


// flatten nested vectors and create Rcpp::vectors for COVER_PROB

// memory for cover probabilities
Rcpp::List RcppListCoverOut; // this is a list of vectors
// 1st element: Rcpp::IntegerVector with fragment IDs as was passed from the R side
// 2nd element: Rcpp::IntegerVector with positions within fragments
// 3rd, 4th and so on: Rcpp::NumericVector with starting probabilities for ftp1, ftp2 and so on
size_t cover_total_size = 0;
for (const auto& v : coverOutFragIDs)
	cover_total_size += v.size();
Rcpp::IntegerVector RcppCoverOutFragIDs(cover_total_size);
offset = 0;
for (const auto& v : coverOutFragIDs) {
	std::copy(v.begin(), v.end(), RcppCoverOutFragIDs.begin() + offset);
	offset += v.size();
}
RcppListCoverOut.push_back(RcppCoverOutFragIDs,"seq");

Rcpp::IntegerVector RcppCoverOutFragPos(cover_total_size);
offset = 0;
for (const auto& v : coverOutFragPos) {
	std::copy(v.begin(), v.end(), RcppCoverOutFragPos.begin() + offset);
	offset += v.size();
}
RcppListCoverOut.push_back(RcppCoverOutFragPos,"pos");

// add flattened start probabilities for each footprint
for(int i=0; i<print_names.size(); ++i){
	Rcpp::NumericVector ftpCoverProbs(cover_total_size,NA_REAL);
	RcppListCoverOut.push_back(ftpCoverProbs,print_names[i]);
}
offset = 0;
for(int seq = 0; seq < coverOutProbs.size(); ++seq){
	for(int i = 0; i < print_names.size(); ++i){
		Rcpp::NumericVector ftpProbVec = RcppListCoverOut[i + 2]; // 0 - fragID, 1 - fragPos, 2 - ftp1, 3 - ftp2 etc.
		std::copy(coverOutProbs[seq][i].begin(), coverOutProbs[seq][i].end(), ftpProbVec.begin() + offset);
	}
	offset += coverOutProbs[seq][0].size();
}

// pack output data into Rcpp:List
Rcpp::List output_data;
output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = RcppListStartOut,
                                  Rcpp::Named("COVER_PROB") = RcppListCoverOut);
return(output_data);
}

// 
// /* // [[Rcpp::depends(RcppProgress)]]
//  * */
// bool Predict::Run(int ncpu)
// {
// 	extern bool _VERBOSE_;
// 	
// 	int numberofobjects = BINDING_OBJECTS.Size();
// 	int seq = 0;
// 	
// 	// initial value for partition sums.
// 	// when footprint priors are normalized, i.e. sum of all priors is 1, then initial values for partition sums is always 1.
// 	// we normalize the priors, therefore we set value to 1.
// 	double part_init = 1;
// 	
// #ifdef _OPENMP
// 	omp_set_nested(true);
// 	omp_set_num_threads(ncpu);
// 	if(_VERBOSE_)
// 		Rcpp::Rcout<<"Running prediction with "<<omp_get_max_threads()<<" cpu."<<endl;
// #endif
// 	
// #pragma omp parallel private(seq)
// {
// 	
// #pragma omp for schedule(dynamic)
// 	for(seq=0;seq < SEQUENCES.Size();++seq){
// 		
// 		
// 		int seqlength = SEQUENCES[seq].Size();
// 		// calculate forward partition summ
// 		F[seq][0] = part_init;
// 		vector<double > pf(numberofobjects,1);
// 		for(int pos=1;pos <= seqlength + BINDING_OBJECTS.maxwmlen; ++pos){
// 			double summ=0;
// 			for(int wm=0;wm < numberofobjects; ++wm){
// 				pf[wm] = 1;
// 				int objlen = BINDING_OBJECTS[wm]->len;
// 				
// 				if(BINDING_OBJECTS[wm]->prior > 0){
// 					pf[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq, pos - objlen);
// 					for(int i = pos - objlen + 1; i <= pos - 1; ++i){
// 						
// 						if(i>=0){
// 							pf[wm] *= F[seq][i];
// 						} else {
// 							pf[wm] *= part_init;
// 						}
// 					}
// 					
// 				}
// 				else{
// 					pf[wm] = 0;
// 				}
// 				summ += pf[wm];
// 			}
// 			
// 			F[seq][pos] = 1/summ;
// 			if(pos<=seqlength){
// 				for(int wm=0;wm<print_indexes.size();++wm){
// 					Prob[seq][wm][pos] = pf[print_indexes[wm]];
// 				}
// 			}
// 		}
// 		
// 		// calculate backward partition summ
// 		vector<double > pb(numberofobjects, 1);
// 		R[seq][seqlength + 1] = part_init;
// 		for(int pos = seqlength; pos>=1;--pos){
// 			double summ = 0;
// 			for(int wm = 0;wm < numberofobjects; ++wm){
// 				pb[wm] = 1;
// 				int objlen = BINDING_OBJECTS[wm]->len;
// 				if(BINDING_OBJECTS[wm]->prior > 0){
// 					pb[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq,pos - 1);
// 					for(int i = pos+1; i <= pos+objlen-1; ++i){
// 						if(i<=seqlength + 1){
// 							pb[wm] *= R[seq][i];
// 						} else{
// 							pb[wm] *= part_init;
// 						}
// 					}
// 					
// 				}
// 				else{
// 					pb[wm] = 0;
// 				}
// 				summ += pb[wm];
// 			}
// 			
// 			R[seq][pos] = 1/summ;
// 		}
// 		// сalculate Z = Fn/Rn
// 		
// 		// calculate initial value for Z based on requirement that total coverage at L must be 1
// 		double zsumm =0;
// 		for(int wm = 0;wm < numberofobjects; ++wm){
// 			int objlen = BINDING_OBJECTS[wm]->len;
// 			double wmsumm = 0;
// 			for(int pos=seqlength; pos<= seqlength + objlen - 1; ++pos){
// 				double prod = 1;
// 				for(int j=pos - objlen + 1; j<=seqlength;++j){
// 					prod *= F[seq][j];
// 				}
// 				wmsumm += pow(part_init,pos - seqlength) * prod;
// 			}
// 			zsumm += BINDING_OBJECTS[wm]->prior * wmsumm;
// 		}
// 		
// 		double z_init = 1/zsumm;
// 		
// 		R[seq][seqlength + 1] = z_init;
// 		for(int pos = seqlength;pos >= 1;--pos){
// 			R[seq][pos] = F[seq][pos] * R[seq][pos + 1]/R[seq][pos];
// 		}
// 		// calculate posteriors
// 		for(int pos = 1;pos <= seqlength;++pos){
// 			
// 			for(int wm=0;wm < print_indexes.size(); ++wm){
// 				int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
// 				if(pos + objlen - 1 <= seqlength)
// 					Prob[seq][wm][pos] = Prob[seq][wm][pos + objlen - 1] * F[seq][pos + objlen -1] * R[seq][pos+objlen];
// 				else
// 					Prob[seq][wm][pos] = 0;
// 			}
// 		}
// 	}
// 	
// }
// 
// return(1);
// }
// 
// 

// 
// double Predict::get_coverage_prob_at_pos_name_index(int seq, int pos, int name_index){ //wm is the index in print_names
// 	
// 	if(Prob.empty()){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
// 	}
// 	
// 	if(seq < 0 || seq>SEQUENCES.Size() - 1){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The sequence index is out of range.\n");
// 	}
// 	
// 	if(pos<1 || pos>SEQUENCES[seq].Size()){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The position is out of range.\n");
// 	}
// 	
// 	if(name_index<0 || name_index > print_names.size() -1){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The name_index is out of range.\n");
// 	}
// 	
// 	string name = print_names[name_index];
// 	double coverage=0;
// 	int objlen = BINDING_OBJECTS[names2indexes[name_index][0]]->len;
// 	for(int i=0;i<names2indicesinprobarray[name_index].size();++i){
// 		for(int p = max(pos - objlen + 1,1); p <= pos; ++p){
// 			coverage += Prob[seq][names2indicesinprobarray[name_index][i]][p];
// 		}
// 	}
// 	
// 	if(coverage<0 || coverage>1+1e-5){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! Incorrect value of coverage probability.\n");
// 	}
// 	
// 	return coverage;
// }
// 
// 
// 
// double Predict::get_coverage_prob_at_pos(int seq, int pos, int wm){ //wm is the index in Prob
// 	
// 	if(Prob.empty()){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
// 	}
// 	
// 	if(seq < 0 || seq>SEQUENCES.Size() - 1){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The sequence index is out of range.\n");
// 	}
// 	
// 	if(pos<1 || pos>SEQUENCES[seq].Size()){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The position is out of range.\n");
// 	}
// 	if(wm<0 || wm > Prob[seq].size() -1){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The wm index is out of range.\n");
// 	}
// 	
// 	double coverage=0;
// 	int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
// 	for(int p = max(pos - objlen + 1,1);p <= pos; ++p)
// 		coverage += Prob[seq][wm][p];
// 	
// 	if(coverage < 0 || coverage > 1+1e-5){
// 		Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! Incorrect value of coverage probability.\n");
// 	}
// 	return coverage;
// 	
// }
// 
// 
// 
// Rcpp::List Predict::getStartProbDF(bool report_prediction_in_flanks){
// 	extern bool _VERBOSE_;
// 	
// 	Rcpp::List out_list;
// 	vector<int32_t > seqnames;
// 	vector<int32_t > positions;
// 	// fill seq and pos
// 	for(int seq=0;seq < SEQUENCES.Size();seq++){
// 		int firstDatPos = SEQUENCES[seq]._firstDatpos;
// 		int lastDatPos = SEQUENCES[seq]._lastDatpos;
// 		
// 		// set start pos for reporting
// 		int spos = report_prediction_in_flanks ? 1 : firstDatPos;
// 		int lpos = lastDatPos;
// 		for(int position = spos; position <= lpos; ++position){
// 			seqnames.push_back(SEQUENCES[seq].Name());
// 			positions.push_back(position - firstDatPos + 1);
// 		}
// 	}
// 	
// 	out_list.push_back(Rcpp::wrap(seqnames),"seq");
// 	out_list.push_back(Rcpp::wrap(positions),"pos");
// 	
// 	// fill start probabilities
// 	for(int i=0;i<print_names.size();++i){
// 		vector<double > stprob;
// 		int seq=0;
// 		
// 		for(seq=0;seq < SEQUENCES.Size();seq++){
// 			int firstDatPos = SEQUENCES[seq]._firstDatpos;
// 			int lastDatPos = SEQUENCES[seq]._lastDatpos;
// 			
// 			int spos = report_prediction_in_flanks ? 1 : firstDatPos;
// 			int lpos = lastDatPos;
// 			
// 			for(int position = spos; position <= lpos; ++position){
// 				double totalprob=0;
// 				for(int j=0;j<names2indicesinprobarray[i].size();++j){
// 					totalprob += Prob[seq][names2indicesinprobarray[i][j]][position];
// 				}
// 				stprob.push_back(totalprob);
// 			}
// 		}
// 		
// 		out_list.push_back(Rcpp::wrap(stprob),print_names[i]);
// 	}
// 	
// 	return out_list;
// }
// 
// 
// 
// Rcpp::List Predict::getCoverProbDF(){
// 	
// 	Rcpp::List out_list;
// 	vector<uint32_t > seqnames;
// 	vector<uint32_t > positions;
// 	// fill seq and pos
// 	for(int seq=0;seq < SEQUENCES.Size();seq++){
// 		int firstDatPos = SEQUENCES[seq]._firstDatpos;
// 		int lastDatPos = SEQUENCES[seq]._lastDatpos;
// 		
// 		for(int position = firstDatPos; position <= lastDatPos; ++position){
// 			seqnames.push_back(SEQUENCES[seq].Name());
// 			positions.push_back(position - firstDatPos + 1);
// 		}
// 	}
// 	
// 	out_list.push_back(Rcpp::wrap(seqnames),"seq");
// 	out_list.push_back(Rcpp::wrap(positions),"pos");
// 	
// 	// fill coverage probabilities
// 	for(int i=0;i<print_names.size();++i){
// 		
// 		vector<double > covprob;
// 		for(int seq=0;seq < SEQUENCES.Size();seq++){
// 			int firstDatPos = SEQUENCES[seq]._firstDatpos;
// 			int lastDatPos = SEQUENCES[seq]._lastDatpos;
// 			
// 			for(int position = firstDatPos; position <= lastDatPos; ++position){
// 				double coverprob=get_coverage_prob_at_pos_name_index(seq,position,i);
// 				covprob.push_back(coverprob);
// 			}
// 		}
// 		out_list.push_back(Rcpp::wrap(covprob),print_names[i]);
// 	}
// 	
// 	return out_list;
// }
// 

// 
// Rcpp::List Predict::getGenomeSummaryDF(){
// 	
// 	// calculate statistics across the whole amplicon (all fragments)
// 	SetGenomeSummary();
// 	Rcpp::List out_list;
// 	const char *fnms[] = {"Prior",
//                        "Expected number of sites",
//                        "Coverage",
//                        "Sites<0.5",
//                        "Sites>=0.5",
//                        "Positions<0.5",
//                        "Positions>=0.5"};
// 	vector<string > field_names(fnms,fnms + 7) ;
// 	out_list.push_back(Rcpp::wrap(field_names),"Statistics");
// 	for(int wm=0;wm<print_names.size();++wm){
// 		out_list.push_back(Rcpp::wrap(genomesummary[wm]),print_names[wm]);
// 	}
// 	
// 	return out_list;
// }


void Predict::clear(){
	
	// SEQUENCES.clear();
	// PARAMS.clear();
	// BINDING_OBJECTS.clear();
	// 
	// F.clear();//swap(vector<vector<double> >());
	// R.clear();//(vector<vector<double> >());
	// Prob.clear();//(vector<vector<vector<double> > >());
	// genomesummary.clear();//(vector<vector<double > > ()); //this vector contains expected prior, expected number of sites and expected coverage for the whole genome and for each TF
	// 
	print_indexes.clear();//swap(vector<int >());	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
	names2indexes.clear();//swap(vector<vector<int > >()); // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
	print_names.clear();//swap(vector<string >()); // this array contain names of the objects that will be printed
	
	names2indicesinprobarray.clear();//swap(vector<vector<int > >()); // this array contains map i - index in names to subarray of indices in Prob array
	
}

