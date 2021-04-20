#ifndef _region_data_
#define _region_data_

#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

#include "fragData-class.hpp"
#include "utils_globvars.hpp"
#include "protectStats-class.hpp"

using namespace std;


// container for storing all fragments overlapping a region of interest
class regionData{
	string regChr;
	int regStart;
	int regEnd;
	
	unordered_map<string, fragData > fDataMap; // map: qname -> fragment
	
public:
	regionData();
	regionData(const string& regchr,
            const int& regstart,
            const int& regend);
	~regionData();
	
	bool addFragContextData(const string& bampref,
                         const string& qname,
                         const int& fragstart,
                         const int& fragend,
                         const vector<int >& rposvec,
                         const vector<bool >& protectvec,
                         const vector<bool >& isonplusstrandvec,
                         const bool& isSecond,// data for same positions on the second read are ignored
                         fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
	);
	
	void print();
	
	
	void addCountsToCoocTable(coocCntTable& cTable);
	
	
	void addCountsToProtectStatsTable(protectStatsTable& pStatsTable);
	
	
	// clip fragments to get rid of partial footprints
	int clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                     const double& max_protect_frac); // keeps fragment if percentage of 1s in a fragment is above this number)
	
	
	
	// get set of qnames whose meth profile is not unique for removal
	set<string > getNonUniqueQnames();
	
	// get set of qnames with failed bisulfite conversion
	set<string > getFailedBisConvQnames(double max_bisC_meth,int min_bisC_size);
	
	// remove fragments from the map. returns new size of the map
	int eraseQnames(const set<string >& qnames);
	
	int size();
	
	
	// convert data in region into a Matrix and return
	Rcpp::IntegerMatrix getDataMatrix(fragMapType whichMap);
	
};



#endif