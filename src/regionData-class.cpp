#include "regionData-class.hpp"

regionData::regionData(){};
regionData::regionData(const string& regchr,
           const int& regstart,
           const int& regend){
  regChr = regchr;
  regStart = regstart;
  regEnd = regend;
};
regionData::~regionData(){};

bool regionData::addFragContextData(const string& bampref,
                        const string& qname,
                        const int& fragstart,
                        const int& fragend,
                        const vector<int >& rposvec,
                        const vector<bool >& protectvec,
                        const vector<bool >& isonplusstrandvec,
                        const bool& isSecond,// data for same positions on the second read are ignored
                        fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
){
  // check if qname exists.
  // no - create new fragment, add data and insert into fDataMap
  // yes - add additional data into corresponding fragment
  
  if(fDataMap.find(qname) == fDataMap.end()){ // if not found a query name then insert a new
    fragData frg;
    frg.addData(bampref,
                fragstart,
                fragend,
                rposvec,
                protectvec,
                isonplusstrandvec,
                isSecond,// data for same positions on the second read are ignored
                whichMap);
    fDataMap.insert(make_pair(qname, frg));
    
  }
  else { // if found then add to existing fragment.
    fDataMap[qname].addData(bampref,
                            fragstart,
                            fragend,
                            rposvec,
                            protectvec,
                            isonplusstrandvec,
                            isSecond,// data for same positions on the second read are ignored
                            whichMap);
  }
  return(1);
};

void regionData::print(){
  Rcpp::Rcout<<"=========== Data in "<<regChr<<":"<<regStart<<"-"<<regEnd<<" ==========="<<endl;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    Rcpp::Rcout<<"==== fragment name: "<<it->first<<" ===="<<endl;
    (it->second).print();
  }
  
  Rcpp::Rcout<<"========================================================================"<<endl;
};


void regionData::addCountsToCoocTable(coocCntTable& cTable){
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    (it->second).addCountsToCoocTable(cTable);
  }
};


void regionData::addCountsToProtectStatsTable(protectStatsTable& pStatsTable){
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    (it->second).addCountsToProtectStatsTable(pStatsTable);
  }
  
}


// clip fragments to get rid of partial footprints
int regionData::clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                    const double& max_protect_frac) // keeps fragment if percentage of 1s in a fragment is above this number)
{
  int nFragsClipped = 0;
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    nFragsClipped += (it->second).clipProtectData(clip_until_nbg,max_protect_frac);
  }
  
  return(nFragsClipped);
};


// get set of qnames whose meth profile is not unique for removal
set<string > regionData::getNonUniqueQnames(){
  set<string > nUqnames;
  unordered_set<string > fragdatenc;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    
    string enc = (it->second).getStringEncodingForUniqueness();
    // for testing
    // (it->second).print();
    // Rcpp::Rcout<<enc<<endl;
    
    if(fragdatenc.find(enc) != fragdatenc.end()){ // encoding exists
      nUqnames.insert(it->first); // add qname to set of non-unique fragments
    } else{ // new encoding
      fragdatenc.insert(enc); 
    }
  }
  return(nUqnames);
  
};

// get set of qnames with failed bisulfite conversion
set<string > regionData::getFailedBisConvQnames(double max_bisC_meth,int min_bisC_size){
  set<string > failed_qnames;
  
  for(unordered_map<string, fragData >::iterator it = fDataMap.begin(); it != fDataMap.end();++it){
    double biscmeth = (it->second).get_bisC_meth(min_bisC_size);
    if(biscmeth > max_bisC_meth)
      failed_qnames.insert(it->first);
  }
  return(failed_qnames);
  
};

// remove fragments from the map. returns new size of the map
int regionData::eraseQnames(const set<string >& qnames){
  for(set<string >::iterator it=qnames.begin(); it != qnames.end(); ++it)
    fDataMap.erase(*it);
  return(fDataMap.size());
};

int regionData::size(){
  return(fDataMap.size());
};


// convert data in region into a Matrix and return
Rcpp::IntegerMatrix regionData::getDataMatrix(fragMapType whichMap){
  
  // create Matrix filled with NA_Integer
  Rcpp::IntegerMatrix matr(fDataMap.size(),regEnd - regStart);
  fill(matr.begin(),matr.end(),NA_INTEGER);
  // set column names
  Rcpp::CharacterVector colnm(regEnd - regStart);
  for(int pos=regStart + 1; pos <= regEnd; ++pos)
    colnm[pos - regStart - 1] = to_string(pos);
  Rcpp::colnames(matr) = colnm;
  // set rownames
  vector<string > rownm;
  for(unordered_map<string, fragData >::iterator it=fDataMap.begin(); it != fDataMap.end(); ++it)
    rownm.push_back(it->first);
  Rcpp::rownames(matr) = Rcpp::wrap(rownm);
  
  if(whichMap != allC){
    for(int rowi=0; rowi < rownm.size();++rowi){
      // get positions in a fragment
      vector<int > pos_vec = fDataMap[rownm[rowi]].getPositions(whichMap);
      for(int posi=0; posi < pos_vec.size(); ++posi){
        if(pos_vec[posi] >= regStart && pos_vec[posi] < regEnd)
          matr(rowi,pos_vec[posi] - regStart) = fDataMap[rownm[rowi]].getDataAt(pos_vec[posi],whichMap);
      }
    }
  } else{
    // if data for all C's is requested then loop across all contexts and fill the matrix
    vector<fragMapType > frag_map_vec = {GCH, WCG, bisC, otherC};
    for(int rowi=0; rowi < rownm.size();++rowi){
      // loop across context map type
      for(int mtpi = 0; mtpi < frag_map_vec.size(); ++mtpi){
        // get positions in a fragment
        vector<int > pos_vec = fDataMap[rownm[rowi]].getPositions(frag_map_vec[mtpi]);
        for(int posi=0; posi < pos_vec.size(); ++posi){
          if(pos_vec[posi] >= regStart && pos_vec[posi] < regEnd){
            // get data at position
            int posdat = fDataMap[rownm[rowi]].getDataAt(pos_vec[posi],frag_map_vec[mtpi]);
            // assign to matrix element, HOWEVER invert for context GCH as we store protection values which are 0-unprotected->methylated and 1-protected->unmethylated
            matr(rowi,pos_vec[posi] - regStart) = frag_map_vec[mtpi] == GCH ? 1-posdat : posdat;
          }
          
        }
      }
    }
  }
  
  return(matr);
};
