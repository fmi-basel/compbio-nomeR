#include "protectStats-class.hpp"

// constructors/destructors
protectStatsTable::protectStatsTable(){};
protectStatsTable::~protectStatsTable(){};

bool protectStatsTable::addCounts(int n_total_gch,
               int n_protect_gch,
               int count){
  if(pStatTable.find(n_total_gch) != pStatTable.end()){ // if row with total number of GCH exists
    
    if(pStatTable[n_total_gch].find(n_protect_gch) != pStatTable[n_total_gch].end()){ // if entry with protected number of gch exists
      pStatTable[n_total_gch][n_protect_gch] += count;
    } else{
      pStatTable[n_total_gch].insert(make_pair(n_protect_gch,count));
    }
    
  } else { // if row with total number of GCH DOESN'T exists
    unordered_map<int, int> prot_count_map;
    prot_count_map.insert(make_pair(n_protect_gch,count));
    pStatTable.insert(make_pair(n_total_gch,prot_count_map));
  }
  
  return(1);
};


vector<vector<int > > protectStatsTable::getStatTable(){
  vector<vector<int > > statTable;
  
  for(unordered_map<int, unordered_map<int, int>>::iterator totit = pStatTable.begin(); totit != pStatTable.end(); ++totit){
    for(unordered_map<int, int>::iterator protit = (totit->second).begin(); protit != (totit->second).end(); ++protit){
      vector<int > row_table;
      row_table.push_back(totit->first); // add total number of gch
      row_table.push_back(protit->first); // add number of protected gch
      row_table.push_back(protit->second); // add number of fragments with this combination
      
      statTable.push_back(row_table); // add row into the output table
    }
  }
  
  return(statTable);
  
};
