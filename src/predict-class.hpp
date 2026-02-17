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

public:

    Predict();
    ~Predict();

    void getCoverPosteriors(const vector<vector<double > >& startProb,
                            const DNAbind_obj_vector& ftpModels,
                            vector<vector<double >>& coverProb
    );


    void getOutputVectors(const vector<vector<double > >& ftpNameProbs,
                          const DNAbind_obj_vector& ftpModels,
                          const fragProtectData& seqData, // current protection data sequence
                          const int& startFrom, // index in seqData to start aggregation
                          const int& endAt, // index in seqData until which to perform aggregation (including)
                          const bool& aggrByGroup, // aggregate by group?
                          vector<int32_t >& outFragIDs,
                          vector<int32_t >& outFragPos,
                          vector<vector<double >>& outProbs // matrix to store probablities, aggregated or not
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

    void getPosteriorViterbiFtpConf(const vector<vector<double >>& coverProb,
                                    const DNAbind_obj_vector& ftpModels,
                                    const fragProtectData& seqData, // current protection data sequence
                                    vector<int32_t >& cPVFragPos,
                                    vector<int32_t >& cPVFtpWidth,
                                    vector<string >& cPVFtpName,
                                    vector<string >& cPVFtpGroup,
                                    vector<double >& cPVFtpProb
    );

    void getPosteriorDecodingFtpConf(const vector<vector<double >>& outCoverProb,
                                     const vector<int32_t >& outCoverFragPos,
                                     const DNAbind_obj_vector& ftpModels,
                                     vector<int32_t >& cPDFragPos,
                                     vector<int32_t >& cPDFtpWidth,
                                     vector<string >& cPDFtpName,
                                     vector<string >& cPDFtpGroup,
                                     vector<double >& cPDFtpProb
    );

    Rcpp::List calcStartCoverProbs(const SMFdataset& smfData,
                                   const DNAbind_obj_vector& ftpModels,
                                   const parameters& params,
                                   ftpConfigAlgo ftpCnfAlg,
                                   bool aggrByGroup,
                                   bool keepStartProb,
                                   int ncpu);



};

#endif
