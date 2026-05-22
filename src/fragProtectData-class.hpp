#ifndef _fragprotectdata_h_
#define _fragprotectdata_h_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <cstdint>
#include <math.h>
#include <Rcpp.h>
#include "utils_globvars.hpp"

using namespace std;


class fragProtectData{

	// all positions within fragments are 1 - based!
public:
	// public data
	uint32_t _fragID; // the index of the fragment that was passed from the R code
	uint32_t _size; // in the previously used Sequence class this was the length of extended fragment
	          // after we padded NAs from left and right sides of length of maximum WM (footprint) size
	          // this was done to be able to calculate probabilities at the edges of molecule
	          // this parameter used in Predict.Run() and other classes.
	          // we have to keep this because memory allocation in Predict class takes this into account.
	uint32_t _firstDatpos; // a position within the extended (by NAs) sequence with first actuall data point
	uint32_t _lastDatpos;  // similarly this is the last position within extended (by NAs) sequence with data points

	vector<double> _modProbVec; // expanded vector with mod prob in [0,1]; NA positions encoded as -1.0

	// public functions
	// constructors/copying/destructors
	fragProtectData();
	fragProtectData(const fragProtectData & s);
	fragProtectData(const uint32_t fragID,
                 const vector<uint32_t>& fragPosVec, // input positions fragPosVec must be 1-based
                 const vector<double>& modProbVec,
                 int maxWMlen);
	~fragProtectData();


	uint32_t Size() const;    // return size of the fragmentData including the padding by NAs by maxWMlen
	uint32_t Name() const; // return fragment ID.
	// int fragLength const;// return actual genomic size of the data, i.e. maximum position with protection value

	const double operator [](uint32_t i) const;

	//vector<uint8_t > subseq(int start,int end);
	fragProtectData & operator = (const fragProtectData & other);
};


#endif
