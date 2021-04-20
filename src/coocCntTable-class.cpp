#include "coocCntTable-class.hpp"


// constructors
coocCntTable::coocCntTable(){};
coocCntTable::coocCntTable(int max_spac){
  max_spacing = max_spac;
  for(int i = 0; i < max_spacing; ++i){
    vector<int > tmp(4,0);
    count_table.push_back(tmp);
  }
};
coocCntTable::~coocCntTable(){};
// getters
int coocCntTable::get_max_spacing(){
  return(max_spacing);
};
vector<vector<int> > coocCntTable::get_count_table(){
  return(count_table);
};

int coocCntTable::get_element(int row, int col){
  return(count_table[row][col]);
}

// add cooc counts to table
bool coocCntTable::addCount(int spacing, // this is distance between positions,
              // where spacing = 0 corresponds to distance 0, i.e. the same position
              int cooc_type_idx, // this is an index of a column in the table,
              // where 0 - N00, 1 - N01, 2 - N10, 3 - N11
              int count){
  if(spacing < max_spacing && spacing >=0){
    count_table[spacing ][cooc_type_idx] += count;
  }
  return(1);
};
