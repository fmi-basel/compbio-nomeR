#include "predict-class.hpp"


Predict::~Predict(){

}

Predict::Predict()
{

}


void Predict::getOutputVectors(const vector<vector<double > >& ftpNameProbs,
                               const DNAbind_obj_vector& ftpModels,
                               const fragProtectData& seqData, // current protection data sequence
                               const int& startFrom, // index in seqData to start aggregation
                               const int& endAt, // index in seqData until which to perform aggregation (including)
                               const bool& aggrByGroup, // aggregate by group?
                               vector<int32_t >& outFragIDs,
                               vector<int32_t >& outFragPos,
                               vector<vector<double >>& outProbs // matrix to store probablities, aggregated or not
){
    int firstDatPos = seqData._firstDatpos;
    int lastDatPos = seqData._lastDatpos;
    //// 1. fill fragIDs and fragPos
    for(int position = startFrom; position <= endAt; ++position){
        outFragIDs.push_back(seqData.Name());
        outFragPos.push_back(position - firstDatPos + 1);
    }

    //// 2. fill probabilities aggregated for each footprint group if aggrByGroup is TRUE
    if(aggrByGroup){
        int nFtpGroups = ftpModels.groups.size();
        for(int igroup = 0; igroup < nFtpGroups; ++igroup){
            // get name for the current ftp group
            string groupName = ftpModels.groups[igroup];
            // get indices for current ftp group
            const vector<int >& groupFtpIndices = ftpModels.getGroupIndexVec(groupName);
            vector<double > cGroupProbs(endAt - startFrom + 1,0);
            for(int position = startFrom; position <= endAt; ++position){
                // summ across probablities associated with current ftp
                double totalprob=0;
                for(int gFtpi = 0; gFtpi < groupFtpIndices.size(); ++gFtpi){
                    totalprob += ftpNameProbs[groupFtpIndices[gFtpi]][position + 1];
                }
                cGroupProbs[position - startFrom] = totalprob;
            }
            outProbs.push_back(cGroupProbs);
        }
    } else{
        for(int ftpIdx = 0; ftpIdx < ftpModels.Size(); ++ftpIdx){
            vector<double > cFtpProbs(endAt - startFrom + 1,0);
            for(int position = startFrom; position <= endAt; ++position){
                cFtpProbs[position - startFrom] = ftpNameProbs[ftpIdx][position + 1];
            }
            outProbs.push_back(cFtpProbs);
        }
    }
}

void Predict::getCoverPosteriors(const vector<vector<double > >& startProb,
                                 const DNAbind_obj_vector& ftpModels,
                                 vector<vector<double >>& coverProb
){
    size_t nFtps = ftpModels.Size();
    size_t problen = startProb[0].size();

    for(int iFtp = 0; iFtp < nFtps; ++iFtp){
        vector<double > cFtpCoverProb(problen,0);
        cFtpCoverProb[0] = startProb[iFtp][0];
        int objlen = ftpModels[iFtp]->len;
        for(int pos = 1; pos < problen; ++pos){
            // add the next start prob to the value at the preceding position
            cFtpCoverProb[pos] = cFtpCoverProb[pos - 1] + startProb[iFtp][pos];
            // subtract start probability at position pos - objlen
            if(pos - objlen >= 0)
                cFtpCoverProb[pos] -= startProb[iFtp][pos - objlen];
        }
        coverProb.push_back(cFtpCoverProb);
    }
}



// Posterior-Viterbi on ftpName coverProbs
void Predict::getPosteriorViterbiFtpConf(const vector<vector<double >>& coverProb,
                                         const DNAbind_obj_vector& ftpModels,
                                         const fragProtectData& seqData, // current protection data sequence
                                         vector<int32_t >& cPVFragPos,
                                         vector<int32_t >& cPVFtpWidth,
                                         vector<string >& cPVFtpName,
                                         vector<string >& cPVFtpGroup,
                                         vector<double >& cPVFtpProb
){
    size_t seqlength = coverProb[0].size(); // length of coverage posteriors vector
    int fDPos = seqData._firstDatpos;
    int lDPos = seqData._lastDatpos;

    // define negative infinity
    // log(0) = -infinity
    const double NEG_INF = -std::numeric_limits<double>::infinity();

    // calculate cumulative logs of coverage posteriors
    // cumSumLogCoverProb[group][pos] is sum of logs in the range [0;pos-1] in the coverProb
    // isValidCumSum[group][pos] is vector keeping validity in case there is any 0 probabilities and to take break cumulative sum
    vector<vector<double >> cumSumLogCoverProb(coverProb.size(), std::vector<double>(seqlength + 1, NEG_INF));
    vector<vector<char >> isValidCumSum(coverProb.size(), std::vector<char>(seqlength + 1, 1));
    for (size_t iFtp = 0; iFtp < coverProb.size(); ++iFtp) {
        double cumsum = 0.0;

        // Position 0 = empty prefix, always valid
        cumSumLogCoverProb[iFtp][0] = 0.0;
        isValidCumSum[iFtp][0] = 1;

        for (size_t pos = 1; pos <= seqlength; ++pos) {
            double p = coverProb[iFtp][pos - 1];
            if (p > 0.0) {
                cumsum += log(p);
                isValidCumSum[iFtp][pos] = 1;
            } else {
                // mark position for iFtp as invalid
                isValidCumSum[iFtp][pos] = 0;
            }

            cumSumLogCoverProb[iFtp][pos] = cumsum;
        }
    }

    vector<double > lFMaxProb(seqlength + 1,0); // vector that keeps maximum configuration log probabilites
    vector<int32_t > ftpEndsTrace(seqlength + 1,-1); // vector containing footprint index with maximum log probability to trace back configuration


    // Viterbi pass using log of cover posteriors as scores
    // pos is index in cumSumLogCoverProb that is pos - 1 in the input coverProb vector
    for(int pos = 1; pos <= seqlength; ++pos){
        double maxLogProb = NEG_INF;
        int bestFtpIdx = -1;

        // go across footprints to find maximum
        for(size_t iFtp = 0; iFtp < cumSumLogCoverProb.size(); ++iFtp){
            int objlen = ftpModels[iFtp]->len; // length of the current footprint
            double curF;
            if(pos - objlen >= 0){
                // cumSumLogCoverProb[group][pos] is sum of logs in the range [0;pos - 1] of the vector ftpGroupCoverProb
                double logScore;
                if(isValidCumSum[iFtp][pos] && isValidCumSum[iFtp][pos - objlen])
                    logScore = cumSumLogCoverProb[iFtp][pos] - cumSumLogCoverProb[iFtp][pos - objlen];
                else
                    logScore = NEG_INF;
                curF = lFMaxProb[pos - objlen] + logScore;
            } else{
                curF = NEG_INF;
            }
            if(curF > maxLogProb){
                maxLogProb = curF;
                bestFtpIdx = iFtp;
            }
        }
        // store best values
        //Rcpp::Rcout<<"Best values: pos="<<posVecCoverProb[pos-1]<<"; maxLogProb="<<maxLogProb<<"; bestFtp="<<ftpModels[bestFtpIdx]->name<<endl;
        lFMaxProb[pos] = maxLogProb;
        ftpEndsTrace[pos] = bestFtpIdx;
    }

    // trace back and construct the best configuration
    int pos = seqlength;
    while(pos + 1 >=  fDPos){
        int bestFtpIndex = ftpEndsTrace[pos];
        int bestFtpLen = ftpModels[bestFtpIndex]->len;
        if(pos - bestFtpLen <= lDPos){ // add if start of footprint is below lDPos
            cPVFragPos.push_back(pos - bestFtpLen - fDPos); // from pos - bestFtpLen + 1 we subtract +1 to get back to indexing of coverProb and shift by fDPos
            cPVFtpWidth.push_back(bestFtpLen);
            cPVFtpName.push_back(ftpModels[bestFtpIndex]->name);
            cPVFtpGroup.push_back(ftpModels[bestFtpIndex]->group);
            // probability that we report for this algorithm is geometric mean of ftp coverages
            //int bestGroupIdx = ftpGroupBestIdx[pos];
            double gMeanCover = exp((cumSumLogCoverProb[bestFtpIndex][pos] - cumSumLogCoverProb[bestFtpIndex][pos - bestFtpLen])/bestFtpLen);
            cPVFtpProb.push_back(gMeanCover);
        }
        pos = pos - bestFtpLen;
    }
}


// Posterior-Decoding ftp configrations
void Predict::getPosteriorDecodingFtpConf(const vector<vector<double >>& outCoverProb,
                                          const vector<int32_t >& outCoverFragPos,
                                          const DNAbind_obj_vector& ftpModels,
                                          vector<int32_t >& cPDFragPos,
                                          vector<int32_t >& cPDFtpWidth,
                                          vector<string >& cPDFtpName,
                                          vector<string >& cPDFtpGroup,
                                          vector<double >& cPDFtpProb
){
    size_t seqlength = outCoverProb[0].size(); // length of coverage posteriors vector

    // define whether the outCoverProb are aggregated by group
    size_t nFtpGroups = ftpModels.getGroupsSize();
    size_t nFtpNames = ftpModels.Size();

    vector<string > ftpNames;
    vector<string > ftpGroupNames;

    if(outCoverProb.size() == nFtpGroups){ // probabilities are aggregated by group
        ftpGroupNames = ftpModels.groups;
        ftpNames = ftpModels.groups; // we register ftpGroups in the column for ftpNames as well
    } else if(outCoverProb.size() == nFtpNames){ // probabilities are NOT aggregated by ftpGroup
        for(size_t iFtp = 0; iFtp < nFtpNames; ++iFtp){
            ftpNames.push_back(ftpModels[iFtp]->name);
            ftpGroupNames.push_back(ftpModels[iFtp]->group);
        }
    } else {
        Rcpp::stop("ERROR:getPosteriorDecodingFtpConf: Number of rows in outCoverProb must equal to number of footprint groups or names.");
    }


    // define negative infinity
    // log(0) = -infinity
    const double NEG_INF = -std::numeric_limits<double>::infinity();

    // find ftp with highest posterior coverage for position one
    int curFtpIdx = -1;
    int32_t curFtpStart = outCoverFragPos[0];
    double curFtpCumLogScore = NEG_INF;
    double maxProb = 0;
    for (size_t iFtp = 0; iFtp < outCoverProb.size(); ++iFtp) {
        if(outCoverProb[iFtp][0] > maxProb){
            maxProb = outCoverProb[iFtp][0];
            curFtpIdx = iFtp;
            curFtpCumLogScore = log(outCoverProb[iFtp][0]);
        }
    }

    // for each position get ftp with maximum posterior coverage and get intervals of continous ftp coverage
    for(int posI = 1; posI < seqlength; ++posI){
        maxProb = 0;
        int maxFtpIdx = -1;

        // find a footprint with maximum coverage posterior
        for (size_t iFtp = 0; iFtp < outCoverProb.size(); ++iFtp) {
            if(outCoverProb[iFtp][posI] > maxProb){
                maxProb = outCoverProb[iFtp][posI];
                maxFtpIdx = iFtp;
            }
        }

        // continue, if the ftp name is the same as at previous position
        if(maxFtpIdx == curFtpIdx){
            curFtpCumLogScore += log(maxProb);
        } else{  // register the segment and start a new segment
            // record positions of the segment
            cPDFragPos.push_back(curFtpStart);
            int32_t ftpWidth = outCoverFragPos[posI - 1] - curFtpStart + 1;
            cPDFtpWidth.push_back(ftpWidth);

            cPDFtpName.push_back(ftpNames[curFtpIdx]);
            cPDFtpGroup.push_back(ftpGroupNames[curFtpIdx]);
            // probability that we report for this algorithm is the geometric mean of ftp coverages
            double gMeanCover = exp(curFtpCumLogScore/ftpWidth);
            cPDFtpProb.push_back(gMeanCover);

            // start a new segment
            curFtpIdx = maxFtpIdx;
            curFtpStart = outCoverFragPos[posI];
            curFtpCumLogScore = log(maxProb);
        }
    }
}


// Viterbi alogirthm to get configuration of footprints with maximum posterior probability
void Predict::getViterbiMAPftpConf(const vector<vector<double > >& ftpModelsScores,
                                   const vector<vector<double > >& startProb,
                                   const DNAbind_obj_vector& ftpModels,
                                   const int& fDPos, // firstDatPos
                                   const int& lDPos, // lastDatPos
                                   vector<int32_t >& cVitFragPos,
                                   vector<int32_t >& cVitFtpWidth,
                                   vector<string >& cVitFtpName,
                                   vector<string >& cVitFtpGroup,
                                   vector<double >& cVitFtpProb){

    //size_t probVecLen = startProb[0].size(); //length of the probability vector
    size_t seqlength = ftpModelsScores[0].size(); // actuall length of extended sequence
    size_t nFtps = ftpModels.Size(); // number of footprints
    // define negative infinity
    // log(0) = -infinity
    const double NEG_INF = -std::numeric_limits<double>::infinity();

    //lastDatPos - firstDatPos + 1
    vector<double > lFMaxProb(seqlength + 1,0); // vector that keeps maximum configuration log probabilites
    vector<int32_t > ftpEndsTrace(seqlength + 1,-1); // vector containing footprint index with maximum log probability to trace back configuration


    // take logs of model scores
    vector<vector<double >> logFtpModelScores(nFtps, std::vector<double>(seqlength, NEG_INF));
    for (size_t w = 0; w < nFtps; ++w) {
        for (size_t i = 0; i < seqlength; ++i) {
            if (ftpModelsScores[w][i] > 0.0)
                logFtpModelScores[w][i] = log(ftpModelsScores[w][i]);
            else
                logFtpModelScores[w][i] = NEG_INF;
        }
    }


    for(int pos = fDPos; pos <= seqlength; ++pos){
        double maxLogProb = NEG_INF;
        int bestFtpIdx = -1;
        // find footprint that maximizes lFMaxProb[pos - len_w] + logP[pos - len_w + 1]
        for(int wm = 0; wm < nFtps; ++wm){
            int objlen = ftpModels[wm]->len;
            double curF = NEG_INF;
            if(pos - objlen >= 0){
                curF = lFMaxProb[pos - objlen] + logFtpModelScores[wm][pos - objlen + 1];
            } else {
                curF = logFtpModelScores[wm][pos - objlen + 1];
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
    while(pos >= fDPos){
        int bestFtpLen = ftpModels[ftpEndsTrace[pos]]->len;
        if(pos - bestFtpLen + 1 <= lDPos){
            cVitFragPos.push_back(pos - bestFtpLen + 1 - fDPos + 1); //  shift by firstDatPos and make 1-based
            cVitFtpWidth.push_back(bestFtpLen);
            cVitFtpName.push_back(ftpModels[ftpEndsTrace[pos]]->name);
            cVitFtpGroup.push_back(ftpModels[ftpEndsTrace[pos]]->group);
            cVitFtpProb.push_back(startProb[ftpEndsTrace[pos]][pos - bestFtpLen + 2]);
        }
        pos = pos - bestFtpLen;
    }

}




// method that calculates start and cover probabilities and returns a Rcpp::List with calculated data.
Rcpp::List Predict::calcStartCoverProbs(const SMFdataset& smfData,
                                        const DNAbind_obj_vector& ftpModels,
                                        const parameters& params,
                                        ftpConfigAlgo ftpCnfAlg,
                                        bool aggrByGroup,
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

    // allocate native C++ vectors for output probabilities
    vector<int32_t > tmp_vec;
    vector<vector<int32_t >> startOutFragIDs(smfData.Size(),tmp_vec); // vector of vectors with fragment IDs. one per seq
    vector<vector<int32_t >> startOutFragPos(smfData.Size(),tmp_vec); // positions within fragments

    vector<vector<int32_t >> coverOutFragIDs(smfData.Size(),tmp_vec); // vector with fragment IDs as was passed from the R side
    vector<vector<int32_t >> coverOutFragPos(smfData.Size(),tmp_vec); // positions within fragments

    vector<vector<vector<double >>> startOutProbs; // vectors of size nFtpGroups, i.e. for each group . per each seq
    vector<vector<vector<double >>> coverOutProbs;

    // allocate vectors for footprint configurations
    vector<vector<int32_t >> ftpConfOutFragIDs(smfData.Size(),tmp_vec);
    vector<vector<int32_t >> ftpConfOutFragPos(smfData.Size(),tmp_vec);
    vector<vector<int32_t >> ftpConfOutFtpWidth(smfData.Size(),tmp_vec);

    vector<string > tmp_str;
    vector<vector<string >> ftpConfOutFtpName(smfData.Size(),tmp_str);
    vector<vector<string >> ftpConfOutFtpGroup(smfData.Size(),tmp_str);
    vector<double > tmp_dbl;
    vector<vector<double >> ftpConfOutFtpProb(smfData.Size(),tmp_dbl);

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
        Rcpp::Rcout<<"Calculating posterior probabilities using "<<omp_get_max_threads()<<" threads."<<endl;
#endif

    Progress prgbar(smfData.Size(), true);
#pragma omp parallel private(seq)
{

#pragma omp for schedule(dynamic)
    for(seq = 0; seq < smfData.Size(); ++seq){
        if (!Progress::check_abort() ) {
            prgbar.increment(); //update progress bar

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
            for(int pos = 0;pos <= seqlength; ++pos){
                for(int wm = 0;wm < nFtpModels; ++wm){
                    int objlen = ftpModels[wm]->len;
                    if(pos + objlen - 1 <= seqlength)
                        Prob[wm][pos] = Prob[wm][pos + objlen - 1] * F[pos + objlen -1] * R[pos+objlen];
                    else
                        Prob[wm][pos] = 0;
                }
            }

            // calculate cover posteriors
            vector<vector<double >> coverProb;
            coverProb.reserve(ftpModels.Size());
            getCoverPosteriors(Prob,
                               ftpModels,
                               coverProb);

            // get output data structure for current sequence for start and cover posteriors
            // startOutFragIDs, startOutFragPos, startOutProbs[ftp]
            // and coverOutFragIDs, coverOutFragPos, coverOutProbs[ftp]
            // for the current molecule
            // NOTE: startOutProbs and coverOutProbs contain aggregated probabilities per group if aggrByGroup is TRUE


            // fill output vectors for START_PROB
            int firstDatPos = smfData[seq]._firstDatpos;
            int lastDatPos = smfData[seq]._lastDatpos;
            // int spos = firstDatPos - maxwmlen + 1; // by default we report calculated posteriors in the left padded region
            // int lpos = lastDatPos;
            int spos = firstDatPos; // ignore padded regions
            int lpos = lastDatPos;
            vector<int32_t > currSeqStartOutFragIDs;
            currSeqStartOutFragIDs.reserve(lpos - spos + 1);
            vector<int32_t > currSeqStartOutFragPos;
            currSeqStartOutFragPos.reserve(lpos - spos + 1);
            vector<vector<double >> currSeqStartOutProbs;
            if(aggrByGroup)
                currSeqStartOutProbs.reserve(ftpModels.groups.size());
            else
                currSeqStartOutProbs.reserve(ftpModels.Size());

            getOutputVectors(Prob,
                             ftpModels,
                             smfData[seq], // current protection data sequence
                                    spos, // index in seqData to start aggregation
                                    lpos, // index in seqData until which to perform aggregation (including)
                                    aggrByGroup, // aggregate by group?
                                    currSeqStartOutFragIDs,
                                    currSeqStartOutFragPos,
                                    currSeqStartOutProbs // matrix to store probablities, aggregated or not
            );


            // fill output vectors for COVER_PROB
            vector<int32_t > currSeqCoverOutFragIDs;
            currSeqCoverOutFragIDs.reserve(lpos - spos + 1);
            vector<int32_t > currSeqCoverOutFragPos;
            currSeqCoverOutFragPos.reserve(lpos - spos + 1);
            vector<vector<double >> currSeqCoverOutProbs;
            if(aggrByGroup)
                currSeqCoverOutProbs.reserve(ftpModels.groups.size());
            else
                currSeqCoverOutProbs.reserve(ftpModels.Size());
            getOutputVectors(coverProb,
                             ftpModels,
                             smfData[seq], // current protection data sequence
                                    spos, // index in seqData to start aggregation
                                    lpos, // index in seqData until which to perform aggregation (including)
                                    aggrByGroup, // aggregate by group?
                                    currSeqCoverOutFragIDs,
                                    currSeqCoverOutFragPos,
                                    currSeqCoverOutProbs // matrix to store probablities, aggregated or not
            );


            // get footprint decoding using chosen algorithm
            vector<int32_t > currFtpConfOutFragPos;
            vector<int32_t > currFtpConfOutFtpWidth;

            vector<string > currFtpConfOutFtpName;
            vector<string > currFtpConfOutFtpGroup;

            vector<double > currFtpConfOutFtpProb;


            switch(ftpCnfAlg) {
            case VITERBI:{

                getViterbiMAPftpConf(ftpModelsScores,
                                     Prob,
                                     ftpModels,
                                     firstDatPos, // TODO: review this! probably it must be changed to spos
                                     lastDatPos,  // and this to lpos
                                     currFtpConfOutFragPos,
                                     currFtpConfOutFtpWidth,
                                     currFtpConfOutFtpName,
                                     currFtpConfOutFtpGroup,
                                     currFtpConfOutFtpProb);
                break;
            }
            case POSTERIORVITERBI:{

                getPosteriorViterbiFtpConf(coverProb,
                                           ftpModels,
                                           smfData[seq],
                                                  currFtpConfOutFragPos,
                                                  currFtpConfOutFtpWidth,
                                                  currFtpConfOutFtpName,
                                                  currFtpConfOutFtpGroup,
                                                  currFtpConfOutFtpProb
                );
                break;
            }
            case POSTERIORDECODING:{

                getPosteriorDecodingFtpConf(currSeqCoverOutProbs,
                                            currSeqCoverOutFragPos,
                                            ftpModels,
                                            currFtpConfOutFragPos,
                                            currFtpConfOutFtpWidth,
                                            currFtpConfOutFtpName,
                                            currFtpConfOutFtpGroup,
                                            currFtpConfOutFtpProb
                );
                break;
            }
            }

            // move all into pre-allocated vectors
            // data for start probabilities
            startOutFragIDs[seq] = move(currSeqStartOutFragIDs);
            startOutFragPos[seq] = move(currSeqStartOutFragPos);
            startOutProbs[seq] = move(currSeqStartOutProbs);
            // data for cover probabilities
            coverOutFragIDs[seq] = move(currSeqCoverOutFragIDs);
            coverOutFragPos[seq] = move(currSeqCoverOutFragPos);
            coverOutProbs[seq] = move(currSeqCoverOutProbs);
            // data for footprint configuration
            vector<int32_t > currFtpConfOutFragIDs(currFtpConfOutFragPos.size(),
                                                   smfData[seq].Name());
            ftpConfOutFragIDs[seq] = move(currFtpConfOutFragIDs);
            ftpConfOutFragPos[seq] = move(currFtpConfOutFragPos);
            ftpConfOutFtpWidth[seq] = move(currFtpConfOutFtpWidth);
            ftpConfOutFtpName[seq] = move(currFtpConfOutFtpName);
            ftpConfOutFtpGroup[seq] = move(currFtpConfOutFtpGroup);
            ftpConfOutFtpProb[seq] = move(currFtpConfOutFtpProb);
        }
    }

} // end of omp parallel


// flatten nested vectors and create Rcpp::vectors for START_PROB
// memory for start probabilities
Rcpp::List RcppListStartOut; // this is a list of vectors
// 1st element: Rcpp::IntegerVector with fragment IDs as was passed from the R side
// 2nd element: Rcpp::IntegerVector with positions within fragments
// 3rd, 4th and so on: Rcpp::NumericVector with starting probabilities for ftp1, ftp2 and so on

// 0. set how many columns base on aggrByGroup
int nElems = 0;
vector<string > elemNames;
if(aggrByGroup){
    nElems = nFtpGroups;
    elemNames = ftpModels.groups;
}
else{
    nElems = nFtpModels;
    for(int iftp=0; iftp < nFtpModels; ++iftp){
        elemNames.push_back(ftpModels[iftp]->name);
    }

}


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

// add flattened start probabilities for each footprint name/group (depend on aggrByGroup bool)
for(int iElem = 0; iElem < nElems; ++iElem){
    Rcpp::NumericVector ftpStartProbs(start_total_size,NA_REAL);
    RcppListStartOut.push_back(ftpStartProbs,elemNames[iElem]);
}
offset = 0;
for(int seq = 0; seq < startOutProbs.size(); ++seq){
    for(int iElem = 0; iElem < nElems; ++iElem){
        Rcpp::NumericVector ftpProbVec = RcppListStartOut[iElem + 2]; // 0 - fragID, 1 - fragPos, 2 - ftp1, 3 - ftp2 etc.
        std::copy(startOutProbs[seq][iElem].begin(), startOutProbs[seq][iElem].end(), ftpProbVec.begin() + offset);
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
for(int iElem = 0; iElem < nElems; ++iElem){
    Rcpp::NumericVector ftpCoverProbs(cover_total_size,NA_REAL);
    RcppListCoverOut.push_back(ftpCoverProbs,elemNames[iElem]);
}

offset = 0;
for(int seq = 0; seq < coverOutProbs.size(); ++seq){
    for(int iElem = 0; iElem < nElems; ++iElem){
        Rcpp::NumericVector ftpProbVec = RcppListCoverOut[iElem + 2]; // 0 - fragID, 1 - fragPos, 2 - ftp1, 3 - ftp2 etc.
        std::copy(coverOutProbs[seq][iElem].begin(), coverOutProbs[seq][iElem].end(), ftpProbVec.begin() + offset);
    }
    offset += coverOutProbs[seq][0].size();
}


// flatten nested vectors and create Rcpp::vectors for footprint configurations
Rcpp::List RcppListFtpConfOut; // this is a list of vectors
// 1st element: Rcpp::IntegerVector with fragment IDs as was passed from the R side
// 2nd element: Rcpp::IntegerVector with starts of footprints within fragments
// 3rd element: Rcpp::IntegerVector with widths of footprints
// 4th element: Rcpp::CharacterVector with footprint names
// 5th element: Rcpp::CharacterVector with footprint groups
// 6th element: Rcpp::NumericVector with probabilities of footprints

size_t ftpconf_total_size = 0;
for (const auto& v : ftpConfOutFragIDs)
    ftpconf_total_size += v.size();

Rcpp::IntegerVector RcppFtpConfOutFragIDs(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFragIDs) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFragIDs.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFragIDs,"seq");

Rcpp::IntegerVector RcppFtpConfOutFragPos(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFragPos) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFragPos.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFragPos,"start");


Rcpp::IntegerVector RcppFtpConfOutFtpWidth(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFtpWidth) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFtpWidth.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFtpWidth,"width");


Rcpp::CharacterVector RcppFtpConfOutFtpName(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFtpName) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFtpName.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFtpName,"ftp_name");


Rcpp::CharacterVector RcppFtpConfOutFtpGroup(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFtpGroup) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFtpGroup.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFtpGroup,"ftp_group");


Rcpp::NumericVector RcppFtpConfOutFtpProb(ftpconf_total_size);
offset = 0;
for (const auto& v : ftpConfOutFtpProb) {
    std::copy(v.begin(), v.end(), RcppFtpConfOutFtpProb.begin() + offset);
    offset += v.size();
}
RcppListFtpConfOut.push_back(RcppFtpConfOutFtpProb,"score");

// pack output data into Rcpp:List
Rcpp::List output_data;
output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = RcppListStartOut,
                                  Rcpp::Named("COVER_PROB") = RcppListCoverOut,
                                  Rcpp::Named("FOOTPRINT_CONF") = RcppListFtpConfOut);
return(output_data);
}



