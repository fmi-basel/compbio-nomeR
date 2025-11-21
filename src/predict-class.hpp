#ifndef _predict_hpp_
#define _predict_hpp_

#include "parameters-class.hpp"
#include "utils_globvars.hpp"
#include "DNAbindobj_vector-class.hpp"
#include "fragProtectData-class.hpp"
#include "SMFdataset-class.hpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <time.h>
#include <Rcpp.h>
#include <limits>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
class Predict
{

	// const parameters& PARAMS; // reference to an object containing parameters
	// const DNAbind_obj_vector& BINDING_OBJECTS; // reference to an object containing vector of footprint models as well as background model
	// const SMFdataset& SEQUENCES; // reference to an object containing SMF data

public:
	// constructors/destructor
	// Predict(const SMFdataset& refSmfData,
	//        const DNAbind_obj_vector& refFtp_models,
	//        const parameters& refParams);
	Predict();

	~Predict();

	void getCoverProbsMatrix(const vector<vector<double > >& startProb,
                          const DNAbind_obj_vector& ftpModels,
                          const int& fDPos, // firstDatPos
                          const int& lDPos, // lastDatPos
                          const size_t& nFtpGroups,
                          vector<vector<double >>& aggrCoverOutProbs
	);

	void getViterbiMAPftpConf(const vector<vector<double > >& ftpModelsScores,
                                    const vector<vector<double > >& startProb,
                                    const DNAbind_obj_vector& ftpModels,
                                    const int& fDPos, // firstDatPos
                                    const int& lDPos, // lastDatPos
                                    vector<int32_t >& cVitFragPos,
                                    vector<int32_t >& cVitFtpWidth,
                                    vector<string >& cVitFtpName,
                                    vector<string >& cVitFtpGroup,
                                    vector<double >& cVitFtpProb);


	void getIntervalScheduleFtpConf(const vector<vector<double >>& ftpGroupStartProb,
                                 const vector<int32_t >& posVecStartProb,
                                 const vector<vector<double > >& ftpNameStartProb,
                                 const DNAbind_obj_vector& ftpModels,
                                 const int& fDPos, // firstDatPos
                                 const int& lDPos, // lastDatPos
                                 vector<int32_t >& cIntSchedFragPos,
                                 vector<int32_t >& cIntSchedFtpWidth,
                                 vector<string >& cIntSchedFtpName,
                                 vector<string >& cIntSchedFtpGroup,
                                 vector<double >& cIntSchedFtpProb);


	Rcpp::List calcStartCoverProbs(const SMFdataset& smfData,
                                const DNAbind_obj_vector& ftpModels,
                                const parameters& params,
                                bool report_prediction_in_flanks,
                                int ncpu);



};



#endif
