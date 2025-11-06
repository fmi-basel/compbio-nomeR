#ifndef _smfdataset_hpp_
#define _smfdataset_hpp_

#include "fragProtectData-class.hpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Rcpp.h>

using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

class SMFdataset{
	int _nmolecs; // number of SMF molecules in the dataset
	int _totallength;
	vector<fragProtectData > _data;

public:
	SMFdataset();
	~SMFdataset();
	SMFdataset(const Rcpp::IntegerVector& fragIDs,
            const Rcpp::IntegerVector& fragPos,
            const Rcpp::IntegerVector& protectVec,
            int maxWMlen);

	bool create(const Rcpp::IntegerVector& fragIDs,
             const Rcpp::IntegerVector& fragPos,
             const Rcpp::IntegerVector& protectVec,
             int maxWMlen);


	void Add(const uint32_t fragID,
          const vector<uint32_t>& fragPosVec, // input positions fragPosVec must be 1-based
          const vector<uint8_t>& protectVec,
          int maxWMlen);
	void Add(fragProtectData & frag);


	const fragProtectData & operator[](int index) const;
	int Size() const;
	int TotalLength() const;
	void clear();

	// function for counting occurrences of 0,0; 0,1 etc at spacing S
	vector<vector<uint64_t > > count_freq_for_spacings(int maxSpacing,
                                                    int ncpu) const;

};

#endif
