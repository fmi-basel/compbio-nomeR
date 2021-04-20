#ifndef _protect_stats_hpp_
#define _protect_stats_hpp_

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
using namespace std;


// container class for protection statistics
class protectStatsTable{
	unordered_map<int, unordered_map<int, int>> pStatTable; // map: total GCH -> protected GCH -> number of fragments
	
public:
	// constructors/destructors
	protectStatsTable();
	~protectStatsTable();
	
	bool addCounts(int n_total_gch,
                int n_protect_gch,
                int count);
	
	
	vector<vector<int > > getStatTable();
	
};

#endif
