#ifndef _cooc_cnttable_
#define _cooc_cnttable_

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>


using namespace std;


// container class for co-occurrence counts
class coocCntTable{
	vector<vector<int> > count_table;
	int max_spacing;
public:
	// constructors
	coocCntTable(){};
	coocCntTable(int max_spac){
		max_spacing = max_spac;
		for(int i = 0; i < max_spacing; ++i){
			vector<int > tmp(4,0);
			count_table.push_back(tmp);
		}
	};
	~coocCntTable(){};
	// getters
	int get_max_spacing(){
		return(max_spacing);
	};
	vector<vector<int> > get_count_table(){
		return(count_table);
	};
	
	int get_element(int row, int col){
		return(count_table[row][col]);
	}
	
	// add cooc counts to table
	bool addCount(int spacing, // this is distance between positions,
               // where spacing = 0 corresponds to distance 0, i.e. the same position
               int cooc_type_idx, // this is an index of a column in the table,
               // where 0 - N00, 1 - N01, 2 - N10, 3 - N11
               int count){
		if(spacing < max_spacing && spacing >=0){
			count_table[spacing ][cooc_type_idx] += count;
		}
		return(1);
	};
	
};



#endif