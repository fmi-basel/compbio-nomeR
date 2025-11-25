#ifndef _predict_hpp_
#define _predict_hpp_

#include "parameters-class.hpp"
#include "utils_globvars.hpp"
#include "DNAbindobj_vector-class.hpp"
#include "fragProtectData-class.hpp"
#include "SMFdataset-class.hpp"
#include "ftpSegment-struct.hpp"
#include "ftpConfig-class.hpp"
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

public:

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


	void getPriorityOrderedFtpConf(const vector<vector<double >>& ftpGroupStartProb,
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
                                ftpConfigAlgo ftpCnfAlg,
                                int ncpu);



};



#endif
