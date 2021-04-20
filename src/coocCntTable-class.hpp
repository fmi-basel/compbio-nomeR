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
	coocCntTable();
	coocCntTable(int max_spac);
	~coocCntTable();
	// getters
	int get_max_spacing();
	vector<vector<int> > get_count_table();
	
	int get_element(int row, int col);
	
	// add cooc counts to table
	bool addCount(int spacing, // this is distance between positions,
               // where spacing = 0 corresponds to distance 0, i.e. the same position
               int cooc_type_idx, // this is an index of a column in the table,
               // where 0 - N00, 1 - N01, 2 - N10, 3 - N11
               int count);
	
};



#endif