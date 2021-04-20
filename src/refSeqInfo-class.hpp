#ifndef _refseq_info_
#define _refseq_info_

#include <Rcpp.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "utils_globvars.hpp"

using namespace std;


// container class for info about tri-nucl context in reference sequence
// it searches all C's and G's in reference sequence and tri-nucleotide context
class refSeqInfo{
	string refSeq;
	int refStart, refEnd;
	
	// map<int, short int> GCH_pos_plus; // stores positions of C's in GCH context on + strand as keys and 1 as values
	// map<int, short int> GCH_pos_minus; // stores positions of G's in GCH context on - strand as keys and -1 as values
	set<int > GCH_pos_plus;
	set<int > GCH_pos_minus;
	
	// map<int, short int> WCG_pos_plus; // stores positions of C's in WCG context on + strand as keys and 1 as values
	// map<int, short int> WCG_pos_minus; // stores positions of G's in WCG context on - strand as keys and -1 as values
	set<int > WCG_pos_plus;
	set<int > WCG_pos_minus;
	
	// map<int, short int> bisC_pos_plus; // stores positions of C's in bisC context on + strand as keys and 1 as values
	// map<int, short int> bisC_pos_minus; // stores positions of G's in bisC context on - strand as keys and -1 as values
	set<int > bisC_pos_plus;
	set<int > bisC_pos_minus;
	
	// map<int, short int> otherC_pos_plus; // stores positions of C's in other contexts on + strand as keys and 1 as values
	// map<int, short int> otherC_pos_minus; // stores positions of G's in other contexts on - strand as keys and -1 as values
	set<int > otherC_pos_plus;
	set<int > otherC_pos_minus;
public:
	// constructors/destructors
	refSeqInfo(const string& seqstring, //sequence of a region
            const int& refstart, // 0 - based start of reference sequence
            const int& refend    // 0-based end of reference sequence
	);
	
	~refSeqInfo();
	
	// print
	void printGCH();
	
	void printWCG();
	
	void printbisC();
	
	void printotherC();
	
	// get context for a 0-based genomic position. 
	// returns 1 - GCH, 2 - WCG, 3 - bisC, 4 - other C. Positive on "+" strand, negative on "-"
	// returns 0 if context is not GCH, WCG or bisC
	short int getContextForRefPos(const int& refpos);
	
	// get positions in range on reference for certain context and on certain strand
	vector<int > getPositionsBetween(int left_rpos, int right_rpos,fragMapType whichMap, bool isOnPlus);
	
	
	// get substring in reference sequence
	string getRefSubseq(int rpos, int len);
	
};


#endif