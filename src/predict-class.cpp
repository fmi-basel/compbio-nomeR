#include "predict-class.hpp"


Predict::~Predict(){

}

Predict::Predict()
{

}


void Predict::getCoverProbsMatrix(const vector<vector<double > >& startProb,
                                  const DNAbind_obj_vector& ftpModels,
                                  const int& fDPos, // firstDatPos
                                  const int& lDPos, // lastDatPos
                                  const size_t& nFtpGroups,
                                  vector<vector<double >>& aggrCoverOutProbs
){

	// for each group
	for(size_t igroup = 0; igroup < nFtpGroups; ++igroup){
		// get name for the current ftp group
		string groupName = ftpModels.groups[igroup];

		// 1. calculate coverage probailities for each footprint using recurrence relationship
		// get footprint indices for current ftp group
		const vector<int >& groupFtpIndices = ftpModels.getGroupIndexVec(groupName);

		for(int gFtpi = 0; gFtpi < groupFtpIndices.size(); ++gFtpi){
			int objlen = ftpModels[groupFtpIndices[gFtpi]]->len;
			vector<double > currFtpCovProb(lDPos - fDPos + 1,0);
			// calculate initial probability at the fDPos (which is firstDatPos)
			for(int p = max(fDPos - objlen + 1,1); p <= fDPos; ++p){
				currFtpCovProb[0] += startProb[groupFtpIndices[gFtpi]][p + 1];
			}
			aggrCoverOutProbs[igroup][0] += currFtpCovProb[0];

			// calculate the rest coverage probabilities by adding and subtracting the probabilities at the next and behind positions
			for(int position = fDPos + 1; position <= lDPos; ++position){
				currFtpCovProb[position - fDPos] =
					currFtpCovProb[position - fDPos - 1] -
					startProb[groupFtpIndices[gFtpi]][position - objlen + 1] + // substract probability at position left behind
					startProb[groupFtpIndices[gFtpi]][position + 1];  // add probability at current position

				aggrCoverOutProbs[igroup][position - fDPos] += currFtpCovProb[position - fDPos];
			}

		}

	}

}

// Viterbi alogirthm to get configuration of footprints with maximum posterior probability
void Predict::getViterbiMAPftpConf(const vector<vector<double > >& startProb,
                          const DNAbind_obj_vector& ftpModels,
                          const int& fDPos, // firstDatPos
                          const int& lDPos, // lastDatPos
                          vector<int32_t >& cVitFragPos,
                          vector<int32_t >& cVitFtpWidth,
                          vector<string >& cVitFtpName,
                          vector<string >& cVitFtpGroup,
                          vector<double >& cVitFtpProb){
	size_t probVecLen = startProb[0].size(); //length of the probability vector
	size_t seqlength = startProb[0].size() - 1; // actuall length of extended sequence
	size_t nFtps = ftpModels.Size(); // number of footprints
	vector<double > lFMaxProb(probVecLen,0); // vector that keeps maximum configuration log probabilites
	vector<int32_t > ftpEndsTrace(probVecLen,-1); // vector containing footprint index with maximum log probability to trace back configuration

	// define negative infinity
	// log(0) = -infinity
	const double NEG_INF = -std::numeric_limits<double>::infinity();
	// Convert start probs to log probabilities
	vector<vector<double >> logP(nFtps, std::vector<double>(probVecLen, NEG_INF));
	for (size_t w = 0; w < nFtps; ++w) {
		for (size_t i = 0; i < probVecLen; ++i) {
			if (startProb[w][i] > 0.0)
				logP[w][i] = log(startProb[w][i]);
		}
	}




	for(int pos = 1; pos <= seqlength; ++pos){
		double maxLogProb = -numeric_limits<double>::infinity();
		int bestFtpIdx = -1;
		// find footprint that maximizes lFMaxProb[pos - len_w] + logP[pos - len_w + 1]
		for(int wm = 0; wm < nFtps; ++wm){
			int objlen = ftpModels[wm]->len;
			double curF = 0;
			if(pos - objlen >= 0){
				curF = lFMaxProb[pos - objlen] + logP[wm][pos - objlen + 1];
			} else {
				curF = log(ftpModels[wm]->nonInformPosterior);
				//curF = NEG_INF;
			}
			if(curF > maxLogProb){
				maxLogProb = curF;
				bestFtpIdx = wm;
			}
		}
		lFMaxProb[pos] = maxLogProb;
		ftpEndsTrace[pos] = bestFtpIdx;

	}

	// trace back and construct the best configuration.
	int pos = seqlength;
	while(pos >= fDPos + 1){
		int bestFtpLen = ftpModels[ftpEndsTrace[pos]]->len;
		cVitFragPos.push_back(pos - bestFtpLen + 1 - fDPos); // also shift by firstDatPos
		cVitFtpWidth.push_back(bestFtpLen);
		cVitFtpName.push_back(ftpModels[ftpEndsTrace[pos]]->name);
		cVitFtpGroup.push_back(ftpModels[ftpEndsTrace[pos]]->group);
		cVitFtpProb.push_back(startProb[ftpEndsTrace[pos]][pos - bestFtpLen + 1]);
		pos = pos - bestFtpLen;
	}

}


// method that calculates start and cover probabilities and returns a Rcpp::List with calculated data.
//
// /* // [[Rcpp::depends(RcppProgress)]]
//  * */
Rcpp::List Predict::calcStartCoverProbs(const SMFdataset& smfData,
                                        const DNAbind_obj_vector& ftpModels,
                                        const parameters& params,
                                        bool report_prediction_in_flanks,
                                        int ncpu){
	extern bool _VERBOSE_;

	size_t nFtpModels = ftpModels.Size(); // number of footprint models including background
	size_t nFtpGroups = ftpModels.getGroupsSize(); // number of groups of footprint models
	// probabilities will be aggregated per group
	int maxwmlen = ftpModels.maxwmlen; // maximum size of footprint model;
	int seq = 0;
	// initial value for partition sums.
	// when footprint priors are normalized, i.e. sum of all priors is 1, then initial values for partition sums is always 1.
	// we normalize the priors, therefore we set value to 1.
	double part_init = 1;

	// allocate natice C++ vectors for output probabilities
	vector<int32_t > tmp_vec;
	vector<vector<int32_t >> startOutFragIDs(smfData.Size(),tmp_vec); // vector of vectors with fragment IDs. one per seq
	vector<vector<int32_t >> startOutFragPos(smfData.Size(),tmp_vec); // positions within fragments

	vector<vector<int32_t >> coverOutFragIDs(smfData.Size(),tmp_vec); // vector with fragment IDs as was passed from the R side
	vector<vector<int32_t >> coverOutFragPos(smfData.Size(),tmp_vec); // positions within fragments

	vector<vector<vector<double >>> startOutProbs; // vectors of size nFtpGroups, i.e. for each group . per each seq
	vector<vector<vector<double >>> coverOutProbs;

	// allocate vectors for maximum aposteriory configurations
	vector<vector<int32_t >> viterbiOutFragIDs(smfData.Size(),tmp_vec);
	vector<vector<int32_t >> viterbiOutFragPos(smfData.Size(),tmp_vec);
	vector<vector<int32_t >> viterbiOutFtpWidth(smfData.Size(),tmp_vec);

	vector<string > tmp_str;
	vector<vector<string >> viterbiOutFtpName(smfData.Size(),tmp_str);
	vector<vector<string >> viterbiOutFtpGroup(smfData.Size(),tmp_str);
	vector<double > tmp_dbl;
	vector<vector<double >> viterbiOutFtpProb(smfData.Size(),tmp_dbl);

	for(seq = 0; seq < smfData.Size(); ++seq){
		vector<vector<double >> tmpst;
		vector<vector<double >> tmpcv;
		for(size_t i=0; i < nFtpGroups; ++i){
			vector<double > startTmp;
			vector<double > coverTmp;
			tmpst.push_back(startTmp);
			tmpcv.push_back(coverTmp);
		}
		startOutProbs.push_back(tmpst);
		coverOutProbs.push_back(tmpcv);
	}





#ifdef _OPENMP
	omp_set_nested(true);
	omp_set_num_threads(ncpu);
	if(_VERBOSE_)
		Rcpp::Rcout<<"Running prediction with "<<omp_get_max_threads()<<" cpu."<<endl;
#endif

#pragma omp parallel private(seq)
{

#pragma omp for schedule(dynamic)
	for(seq = 0; seq < smfData.Size(); ++seq){
		int seqlength = smfData[seq].Size();

		// calculate footprint model scores for the current fragment
		vector<vector<double >> ftpModelsScores = ftpModels.getFtpModelScores(smfData[seq]);


		// allocate memory for:
		// F -  forward parition sum
		// R - backward partition sum
		// Prob - probability of footprint ends at position pos

		vector<double > F(seqlength + maxwmlen + 1,1); // allocate memory for forward parition sum
		vector<double > R(seqlength + 2,1); // allocate memory for backward partition sum. it is shorter than F
		vector<double > probPerFtp(seqlength + 1, 0);
		vector<vector<double > > Prob(nFtpModels, probPerFtp);
		// vector<vector<double > > Prob(print_indexes.size(),probPerFtp);

		// calculate forward partition summ
		F[0] = part_init;
		vector<double > pf(nFtpModels,1);
		for(int pos = 1; pos <= seqlength + ftpModels.maxwmlen; ++pos){
			double summ=0;
			for(int wm = 0; wm < nFtpModels; ++wm){
				pf[wm] = 1;
				int objlen = ftpModels[wm]->len;

				if(ftpModels[wm]->prior > 0){

					if(pos - objlen >= 0 && pos - objlen < seqlength)
						pf[wm] = ftpModelsScores[wm][pos - objlen];
					else
						pf[wm] = ftpModels[wm]->prior;
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
			if(pos <= seqlength){
				for(int wm = 0; wm < nFtpModels; ++wm){
					Prob[wm][pos] = pf[wm];
				}
			}
		}

		// calculate backward partition summ
		vector<double > pb(nFtpModels, 1);
		R[seqlength + 1] = part_init;
		for(int pos = seqlength; pos >= 1; --pos){
			double summ = 0;
			for(int wm = 0; wm < nFtpModels; ++wm){
				pb[wm] = 1;
				int objlen = ftpModels[wm]->len;
				if(ftpModels[wm]->prior > 0){

					if(pos - 1 >= 0 && pos + objlen - 1 < seqlength)
						pb[wm] = ftpModelsScores[wm][pos - 1];
					else
						pb[wm] = ftpModels[wm]->prior;
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
		double zsumm = 0;
		for(int wm = 0; wm < nFtpModels; ++wm){
			int objlen = ftpModels[wm]->len;
			double wmsumm = 0;
			for(int pos = seqlength; pos <= seqlength + objlen - 1; ++pos){
				double prod = 1;
				for(int j=pos - objlen + 1; j<=seqlength;++j){
					prod *= F[j];
				}
				wmsumm += pow(part_init,pos - seqlength) * prod;
			}
			zsumm += ftpModels[wm]->prior * wmsumm;
		}

		double z_init = 1/zsumm;

		R[seqlength + 1] = z_init;
		for(int pos = seqlength; pos >= 1; --pos){
			R[pos] = F[pos] * R[pos + 1]/R[pos];
		}
		// calculate start posteriors
		// prev for(int pos = 1;pos <= seqlength; ++pos){//
		for(int pos = 0;pos <= seqlength; ++pos){
			for(int wm = 0;wm < nFtpModels; ++wm){
				int objlen = ftpModels[wm]->len;
				if(pos + objlen - 1 <= seqlength)
					Prob[wm][pos] = Prob[wm][pos + objlen - 1] * F[pos + objlen -1] * R[pos+objlen];
				else
					Prob[wm][pos] = 0;
			}
		}


		// get output data structure for current sequence for start and cover probabilities
		// startOutFragIDs, startOutFragPos, startOutProbs[ftp]
		// and coverOutFragIDs, coverOutFragPos, coverOutProbs[ftp]
		// for the current molecule
		// NOTE: startOutProbs and coverOutProbs contain aggregated probabilities per group


		// fill output vectors for START_PROB
		int firstDatPos = smfData[seq]._firstDatpos;
		int lastDatPos = smfData[seq]._lastDatpos;
		int spos = report_prediction_in_flanks ? 1 : firstDatPos;
		int lpos = lastDatPos;
		vector<int32_t > currSeqStartOutFragIDs;
		vector<int32_t > currSeqStartOutFragPos;
		vector<double > tmpstart;
		vector<vector<double >> currSeqStartOutProbs(nFtpGroups,tmpstart);
		for(int position = spos; position <= lpos; ++position){
			//// 1. fill fragIDs and fragPos
			currSeqStartOutFragIDs.push_back(smfData[seq].Name());
			currSeqStartOutFragPos.push_back(position - firstDatPos + 1);
			//// 2. fill start probabilities aggregated for each footprint group
			for(int igroup = 0; igroup < nFtpGroups; ++igroup){
				// get name for the current ftp group
				string groupName = ftpModels.groups[igroup];
				// get indices for current ftp group
				const vector<int >& groupFtpIndices = ftpModels.getGroupIndexVec(groupName);
				// summ across probablities associated with current ftp
				double totalprob=0;
				for(int gFtpi = 0; gFtpi < groupFtpIndices.size(); ++gFtpi){
					totalprob += Prob[groupFtpIndices[gFtpi]][position + 1];
				}
				currSeqStartOutProbs[igroup].push_back(totalprob);
			}
		}
		startOutFragIDs[seq] = move(currSeqStartOutFragIDs);
		startOutFragPos[seq] = move(currSeqStartOutFragPos);
		startOutProbs[seq] = move(currSeqStartOutProbs);


		// fill output vectors for COVER_PROB. Perhaps, this can be optimized by adding and subtracting start prob at end and beginning of footprint
		vector<int32_t > currSeqCoverOutFragIDs;
		vector<int32_t > currSeqCoverOutFragPos;

		for(int position = firstDatPos; position <= lastDatPos; ++position){
			//// 1. fill fragIDs and fragPos
			currSeqCoverOutFragIDs.push_back(smfData[seq].Name());
			currSeqCoverOutFragPos.push_back(position - firstDatPos + 1);
		}


		coverOutFragIDs[seq] = move(currSeqCoverOutFragIDs);
		coverOutFragPos[seq] = move(currSeqCoverOutFragPos);


		// allocate vectors for coverage probabilities for each group
		vector<double > tmpcov(lastDatPos - firstDatPos + 1,0);;
		vector<vector<double >> currSeqCoverOutProbs(nFtpGroups,tmpcov); // aggregated probabilities across all footprints per group;

		// vector<double > tmpcov(lastDatPos - firstDatPos + 1,0);
		// coverOutProbs[seq].resize(nFtpGroups,tmpcov);

		getCoverProbsMatrix(Prob,
                      ftpModels,
                      firstDatPos,
                      lastDatPos,
                      nFtpGroups,
                      currSeqCoverOutProbs
		);
		coverOutProbs[seq] = move(currSeqCoverOutProbs);

		// get maximum aposteriori comfiguration of footprints using Viterbi algorithm

		vector<int32_t > currViterbiOutFragPos;
		vector<int32_t > currViterbiOutFtpWidth;

		vector<string > currViterbiOutFtpName;
		vector<string > currViterbiOutFtpGroup;

		vector<double > currViterbiOutFtpProb;
		getViterbiMAPftpConf(Prob,
                       ftpModels,
                       firstDatPos,
                       lastDatPos,
                       currViterbiOutFragPos,
                       currViterbiOutFtpWidth,
                       currViterbiOutFtpName,
                       currViterbiOutFtpGroup,
                       currViterbiOutFtpProb);
		vector<int32_t > currViterbiOutFragIDs(currViterbiOutFragPos.size(),smfData[seq].Name());
		viterbiOutFragIDs[seq] = currViterbiOutFragIDs;
		viterbiOutFragPos[seq] = currViterbiOutFragPos;
		viterbiOutFtpWidth[seq] = currViterbiOutFtpWidth;
		viterbiOutFtpName[seq] = currViterbiOutFtpName;
		viterbiOutFtpGroup[seq] = currViterbiOutFtpGroup;
		viterbiOutFtpProb[seq] = currViterbiOutFtpProb;


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

// add flattened start probabilities for each footprint group
for(int igroup = 0; igroup < nFtpGroups; ++igroup){
	Rcpp::NumericVector ftpStartProbs(start_total_size,NA_REAL);
	RcppListStartOut.push_back(ftpStartProbs,ftpModels.groups[igroup]);
}
offset = 0;
for(int seq = 0; seq < startOutProbs.size(); ++seq){
	for(int igroup = 0; igroup < nFtpGroups; ++igroup){
		Rcpp::NumericVector ftpProbVec = RcppListStartOut[igroup + 2]; // 0 - fragID, 1 - fragPos, 2 - ftpGroup1, 3 - ftpGroup2 etc.
		std::copy(startOutProbs[seq][igroup].begin(), startOutProbs[seq][igroup].end(), ftpProbVec.begin() + offset);
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

// add flattened cover probabilities for each footprint group
for(int igroup = 0; igroup < nFtpGroups; ++igroup){
	Rcpp::NumericVector ftpCoverProbs(cover_total_size,NA_REAL);
	RcppListCoverOut.push_back(ftpCoverProbs,ftpModels.groups[igroup]);
}

offset = 0;
for(int seq = 0; seq < coverOutProbs.size(); ++seq){
	for(int igroup = 0; igroup < nFtpGroups; ++igroup){
		Rcpp::NumericVector ftpProbVec = RcppListCoverOut[igroup + 2]; // 0 - fragID, 1 - fragPos, 2 - ftp1, 3 - ftp2 etc.
		std::copy(coverOutProbs[seq][igroup].begin(), coverOutProbs[seq][igroup].end(), ftpProbVec.begin() + offset);
	}
	offset += coverOutProbs[seq][0].size();
}

// viterbiOutFragIDs[seq] = currViterbiOutFragIDs;
// viterbiOutFragPos[seq] = currViterbiOutFragPos;
// viterbiOutFtpWidth[seq] = currViterbiOutFtpWidth;
// viterbiOutFtpName[seq] = currViterbiOutFtpName;
// viterbiOutFtpGroup[seq] = currViterbiOutFtpGroup;
// viterbiOutFtpProb[seq] = currViterbiOutFtpProb;
// flatten nested vectors and create Rcpp::vectors for Viterbi MAP footprint configurations
Rcpp::List RcppListViterbiOut; // this is a list of vectors
// 1st element: Rcpp::IntegerVector with fragment IDs as was passed from the R side
// 2nd element: Rcpp::IntegerVector with starts of footprints within fragments
// 3rd element: Rcpp::IntegerVector with widths of footprints
// 4th element: Rcpp::CharacterVector with footprint names
// 5th element: Rcpp::CharacterVector with footprint groups
// 6th element: Rcpp::NumericVector with start probabilities of footprints

size_t viterbi_total_size = 0;
for (const auto& v : viterbiOutFragIDs)
	viterbi_total_size += v.size();

Rcpp::IntegerVector RcppViterbiOutFragIDs(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFragIDs) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFragIDs.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFragIDs,"seq");

Rcpp::IntegerVector RcppViterbiOutFragPos(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFragPos) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFragPos.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFragPos,"start");

Rcpp::IntegerVector RcppViterbiOutFtpWidth(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFtpWidth) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFtpWidth.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFtpWidth,"width");
// viterbiOutFtpName[seq] = currViterbiOutFtpName;
// viterbiOutFtpGroup[seq] = currViterbiOutFtpGroup;
// viterbiOutFtpProb[seq] = currViterbiOutFtpProb;
Rcpp::CharacterVector RcppViterbiOutFtpName(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFtpName) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFtpName.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFtpName,"ftp_name");

Rcpp::CharacterVector RcppViterbiOutFtpGroup(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFtpGroup) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFtpGroup.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFtpGroup,"ftp_group");

Rcpp::NumericVector RcppViterbiOutFtpProb(viterbi_total_size);
offset = 0;
for (const auto& v : viterbiOutFtpProb) {
	std::copy(v.begin(), v.end(), RcppViterbiOutFtpProb.begin() + offset);
	offset += v.size();
}
RcppListViterbiOut.push_back(RcppViterbiOutFtpProb,"start_prob");



// pack output data into Rcpp:List
Rcpp::List output_data;
output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = RcppListStartOut,
                                  Rcpp::Named("COVER_PROB") = RcppListCoverOut,
                                  Rcpp::Named("VITERBI_CONF") = RcppListViterbiOut);
return(output_data);
}



