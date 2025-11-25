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
#include "ftpSegment-struct.hpp"
#include "predict-class.hpp"

using namespace std;

// Binary Indexed Tree (Fenwick tree) implementation for segment occupancy vector
// NOTE: indexing in the original arr is from [1..n]
struct ftpOccupancyVectorBIT {
	int n;
	vector<int> bit1, bit2;

	ftpOccupancyVectorBIT(int n) : n(n), bit1(n+1, 0), bit2(n+1, 0) {}

	void add(vector<int>& bit, int i, int v) {
		for (; i <= n; i += i & -i) bit[i] += v;
	}

	// add v to range [l, r]
	void range_add(int l, int r, int v) {
		if (l > r) return;
		add(bit1, l, v);
		add(bit1, r+1, -v);
		add(bit2, l, v*(l-1));
		add(bit2, r+1, -v*r);
	}

	int sum(const vector<int>& bit, int i) const {
		int s = 0;
		for (; i > 0; i -= i & -i) s += bit[i];
		return s;
	}

	// prefix sum from 1..i
	int prefix_sum(int i) const {
		return sum(bit1, i) * i - sum(bit2, i);
	}

	// range sum from l..r
	int range_sum(int l, int r) const {
		if (l > r) return 0;
		return prefix_sum(r) - prefix_sum(l-1);
	}
};





// class intended to add footprints and store footprint configuration on a molecule without overlaps
class ftpConfig
{
	vector<int32_t > ftpStartPosVec;    // starting positions of footprints in a configuration
	vector<int32_t > ftpWidthVec;       // widths of footprints
	vector<double > ftpGroupStartProbVec;  // start probability for a footprint group, i.e. aggregate across all footprints in the same group
	vector<double > ftpNameStartProbVec;   // start probability for the particular footprint with width W
	vector<string > ftpNameVec;         // footprint names
	vector<string > ftpGroupVec;        // footprint groups

	ftpOccupancyVectorBIT ftpOccupancyVec;          // occupancy in the current configuration
	                                                // NOTE: indexing in the occupancy vector is [1..n]

	size_t length;
	int totalOccupied;
	bool isFool;


	//int64_t findNextOccupPos(const int64_t startPos);   // find next occupied position starting from _pos

// 	void markOccupPositions(const int64_t start,
//                          const int32_t width);

public:
	ftpConfig(size_t molLength);

	~ftpConfig();

	bool canFitFtp(const ftpSegment& cFtpSegm);

	bool addFtp(const ftpSegment& cFtpSegm);

	void fillConfigVector(vector<int32_t >& cIntSchedFragPos,
                       vector<int32_t >& cIntSchedFtpWidth,
                       vector<string >& cIntSchedFtpName,
                       vector<string >& cIntSchedFtpGroup,
                       vector<double >& cIntSchedFtpProb);

	bool isMoleculeFool(){return isFool;};

};

#endif
