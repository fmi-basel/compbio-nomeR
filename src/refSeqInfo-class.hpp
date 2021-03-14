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
	){
		refSeq = stringToUpper(seqstring);
		refStart = refstart;
		refEnd = refend;
		
		// find all C' and G's and store positions in corresponding map
		for(int pos = 0; pos <= refSeq.length() - 3; ++pos){
			string context = refSeq.substr(pos,3);
			if(GCH_CONTEXTS.find(context) != GCH_CONTEXTS.end()){
				if(GCH_CONTEXTS[context] > 0)
					GCH_pos_plus.insert(refStart + pos + 1);
				else
					GCH_pos_minus.insert(refStart + pos + 1);
			} else if(WCG_CONTEXTS.find(context) != WCG_CONTEXTS.end()){
				if(WCG_CONTEXTS[context] > 0)
					WCG_pos_plus.insert(refStart + pos + 1);
				else
					WCG_pos_minus.insert(refStart + pos + 1);
			} else if(BISC_CONTEXTS.find(context) != BISC_CONTEXTS.end()){
				if(BISC_CONTEXTS[context] > 0)
					bisC_pos_plus.insert(refStart + pos + 1);
				else
					bisC_pos_minus.insert(refStart + pos + 1);
			} else if(context[1] == 'C'){
				otherC_pos_plus.insert(refStart + pos + 1);
			} else if(context[1] == 'G'){
				otherC_pos_minus.insert(refStart + pos + 1);
			}
		}
	};
	
	~refSeqInfo(){};
	
	// print
	void printGCH(){
		set<int>::iterator it = GCH_pos_plus.begin();
		Rcpp::Rcout<<"Positive on '+' strand, negative on '-'"<<endl;
		while(it != GCH_pos_plus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: + "<<endl;
			it++;
		}
		it = GCH_pos_minus.begin();
		while(it != GCH_pos_minus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: - "<<endl;
			it++;
		}
		
	};
	
	void printWCG(){
		set<int>::iterator it = WCG_pos_plus.begin();
		Rcpp::Rcout<<"Positive on '+' strand, negative on '-'"<<endl;
		while(it != WCG_pos_plus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: + "<<endl;
			it++;
		}
		it = WCG_pos_minus.begin();
		while(it != WCG_pos_minus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: - "<<endl;
			it++;
		}
	};
	
	void printbisC(){
		set<int >::iterator it = bisC_pos_plus.begin();
		Rcpp::Rcout<<"Positive on '+' strand, negative on '-'"<<endl;
		while(it != bisC_pos_plus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: + "<<endl;
			it++;
		}
		it = bisC_pos_minus.begin();
		while(it != bisC_pos_minus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: - "<<endl;
			it++;
		}
	};
	
	void printotherC(){
		set<int >::iterator it = otherC_pos_plus.begin();
		Rcpp::Rcout<<"Positive on '+' strand, negative on '-'"<<endl;
		while(it != otherC_pos_plus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: + "<<endl;
			it++;
		}
		it = otherC_pos_minus.begin();
		while(it != otherC_pos_minus.end())
		{
			Rcpp::Rcout<<"refpos: "<<*it<<"; strand: - "<<endl;
			it++;
		}
		
	};
	
	// get context for a 0-based genomic position. 
	// returns 1 - GCH, 2 - WCG, 3 - bisC, 4 - other C. Positive on "+" strand, negative on "-"
	// returns 0 if context is not GCH, WCG or bisC
	short int getContextForRefPos(const int& refpos){
		if(GCH_pos_plus.find(refpos) != GCH_pos_plus.end())
			return(1);
		else if(GCH_pos_minus.find(refpos) != GCH_pos_minus.end())
			return(-1);
		else if(WCG_pos_plus.find(refpos) != WCG_pos_plus.end())
			return(2);
		else if(WCG_pos_minus.find(refpos) != WCG_pos_minus.end())
			return(-2);
		else if(bisC_pos_plus.find(refpos) != bisC_pos_plus.end())
			return(3);
		else if(bisC_pos_minus.find(refpos) != bisC_pos_minus.end())
			return(-3);
		else if(otherC_pos_plus.find(refpos) != otherC_pos_plus.end())
			return(4);
		else if(otherC_pos_minus.find(refpos) != otherC_pos_minus.end())
			return(-4);
		else
			return(0);
	};
	
	// get positions in range on reference for certain context and on certain strand
	vector<int > getPositionsBetween(int left_rpos, int right_rpos,fragMapType whichMap, bool isOnPlus){
		
		// select which map to get the positions from
		set<int > *mdat; // pointer to methylation/protection map
		switch(whichMap) {
		case GCH:
		{
			if(isOnPlus)
				mdat = &GCH_pos_plus;
			else
				mdat = &GCH_pos_minus;
			break;
		}
		case WCG:
		{
			if(isOnPlus)
				mdat = &WCG_pos_plus;
			else
				mdat = &WCG_pos_minus;
			break;
		}
		case bisC:
		{
			if(isOnPlus)
				mdat = &bisC_pos_plus;
			else
				mdat = &bisC_pos_minus;
			break;
		}
		case otherC:
			if(isOnPlus)
				mdat = &otherC_pos_plus;
			else
				mdat = &otherC_pos_minus;
			break;
		}
		
		if(left_rpos < refStart || right_rpos > refEnd){
			Rcpp::Rcout<<"refSeqInfo::getGCHpositionsBetween: WARNING: desired position left_rpos="<<left_rpos<<" or right_rpos="<<right_rpos<<" is outside of reference sequence ["<<refStart<<";"<<refEnd<<"]. Ignoring everything outside of this range."<<endl;
		}
		
		set<int >::iterator lit = mdat->lower_bound(left_rpos);
		set<int >::iterator rit = mdat->upper_bound(right_rpos + 1);
		vector<int > rpos_in_range;
		for(set<int >::iterator it = lit; it != rit; ++it){
			rpos_in_range.push_back(*it);
		}
		return(rpos_in_range);
	};
	
	
	// get substring in reference sequence
	string getRefSubseq(int rpos, int len){
		if(rpos < refStart || rpos + len > refEnd){
			
			Rcpp::Rcout<<"refSeqInfo::getRefSubseq: WARNING: requested position and length are outside of reference sequence: rpos="<<rpos<<"; len="<<len<<"; refSeq range ["<<refStart<<";"<<refEnd<<"]. Returning N."<<endl;
			return(string("N"));
		}
		return(refSeq.substr(rpos - refStart,len));
	};
	
};


#endif