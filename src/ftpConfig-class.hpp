#ifndef _ftpconfig_hpp_
#define _ftpconfig_hpp_

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

using namespace std;

// class intended to add footprints and store footprint configuration on a molecule without overlaps
class ftpConfig
{
	vector<int64_t > ftpStartPosVec;    // starting positions of footprints in a configuration
	vector<int32_t > ftpWidthVec;       // widths of footprints
	vector<double > ftpGroupStartProbVec;  // start probability for a footprint group, i.e. aggregate across all footprints in the same group
	vector<double > ftpNameStartProbVec;   // start probability for the particular footprint with width W
	vector<string > ftpNameVec;         // footprint names
	vector<string > ftpGroupVec;        // footprint groups

	vector<uint8_t> _occupancy;          // occupancy in the current configuration


	int64_t findNextOccupPos(const int64_t startPos);   // find next occupied position starting from _pos

	void markOccupPositions(const int64_t start,
                         const int32_t width);




public:
	Predict(size_t molLength);

	~Predict();

	bool addFtp(const int64_t& _ftpStart,
             const int32_t& _ftpWidth,
             const double& _ftpGroupStartProb,
             const double& _ftpNameStartProb,
             const string& _ftpName,
             const string& _ftpGroup);

};

#endif
