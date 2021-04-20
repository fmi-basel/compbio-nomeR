#include "fragData-class.hpp"


fragData::fragData(){
  fragStart = -1;
  fragEnd = -1;
  prefix = "";
};

// constructors/destructors
fragData::fragData(const string& bampref,
         const int& fragstart,
         const int& fragend,
         const vector<int >& gch_pos,
         const vector<bool >& gch_protect,
         const vector<bool >& gch_onplusstrand,
         const vector<int >& wcg_pos,
         const vector<bool >& wcg_protect,
         const vector<bool >& wcg_onplusstrand,
         const vector<int >& bisc_pos,
         const vector<bool >& bisc_protect,
         const vector<bool >& bisc_onplusstrand,
         const vector<int >& otherc_pos,
         const vector<bool >& otherc_protect,
         const vector<bool >& otherc_onplusstrand
){
  
  // set coordinates
  fragStart = fragstart;
  fragEnd = fragend;
  prefix = bampref;
  // add GCH data
  addData(bampref,fragstart, fragend, gch_pos, gch_protect, gch_onplusstrand, 0,GCH);
  // add WCG data
  addData(bampref,fragstart, fragend, wcg_pos, wcg_protect, wcg_onplusstrand, 0,WCG);
  // add bisC data
  addData(bampref,fragstart, fragend, bisc_pos, bisc_protect, bisc_onplusstrand, 0,bisC);
  // add otherC data
  addData(bampref,fragstart, fragend, otherc_pos, otherc_protect, otherc_onplusstrand, 0,otherC);
};

fragData::~fragData(){};

// get vector of positions in a map
vector<int > fragData::getPositions(fragMapType whichMap){
  // select which map to get the data from
  map<int , bool > *mdat; // pointer to methylation/protection map
  
  switch(whichMap) {
  case GCH:
    mdat = &gchProtect;
    break;
  case WCG:
    mdat = &wcgMeth;
    break;
  case bisC:
    mdat = &bisCMeth;
    break;
  case otherC:
    mdat = &otherCMeth;
    break;
  }
  vector<int > pos_vec;
  
  for(map<int, bool>::iterator it = mdat->begin(); it != mdat->end(); ++it){
    pos_vec.push_back(it->first);
  }
  return(pos_vec);
  
}


// get data at position
bool fragData::getDataAt(int rpos,
               fragMapType whichMap){
  
  // select which map to get the data from
  map<int , bool > *mdat; // pointer to methylation/protection map
  
  switch(whichMap) {
  case GCH:
    mdat = &gchProtect;
    break;
  case WCG:
    mdat = &wcgMeth;
    break;
  case bisC:
    mdat = &bisCMeth;
    break;
  case otherC:
    mdat = &otherCMeth;
    break;
  }
  if(mdat->find(rpos) == mdat->end()){
    Rcpp::Rcout<<"Key "<<rpos<<" was not found"<<endl;
    Rcpp::stop("");
  }
  return((*mdat)[rpos]);
};

// get isOnPlus
bool fragData::isDataOnPlusStrand(int rpos,
                        fragMapType whichMap){
  
  // select which map to get the data from
  map<int , bool > *sdat; // pointer to strand map
  switch(whichMap) {
  case GCH:
    sdat = &gchIsOnPlusStrand;
    break;
  case WCG:
    sdat = &wcgIsOnPlusStrand;
    break;
  case bisC:
    sdat = &bisCIsOnPlusStrand;
    break;
  case otherC:
    sdat = &otherCIsOnPlusStrand;
    break;
  }
  
  if(sdat->find(rpos) == sdat->end()){
    Rcpp::Rcout<<"Key "<<rpos<<" was not found"<<endl;
    Rcpp::stop("");
  }
  return((*sdat)[rpos]);
};


// add data into the map 
bool fragData::addData(const string& bampref,
             const int& fragstart,
             const int& fragend,
             const vector<int >& rposvec,
             const vector<bool >& protectvec,
             const vector<bool >& isonplusstrandvec,
             const bool& isSecond,// data for same positions on the second read are ignored
             fragMapType whichMap // must be enumerator GCH (protect), WCG (endogenous), bisC or otherC
){ 
  if(rposvec.size() != protectvec.size()|| rposvec.size() != isonplusstrandvec.size()){
    Rcpp::Rcout<<"fragProtectData::addData: sizes of rposvec, protectvec and isonplusstrandvec do not match.";
    return(0);
  }
  
  // change fragment coordinates if neccessary
  if(fragStart < 0 || fragstart < fragStart)
    fragStart = fragstart;
  if(fragEnd < 0 || fragend > fragEnd)
    fragEnd = fragend;
  prefix = bampref;
  
  // select which map to add the data to
  map<int , bool > *mdat; // pointer to methylation/protection map
  map<int , bool > *sdat; // pointer to strand map
  switch(whichMap) {
  case GCH:
    mdat = &gchProtect;
    sdat = &gchIsOnPlusStrand;
    break;
  case WCG:
    mdat = &wcgMeth;
    sdat = &wcgIsOnPlusStrand;
    break;
  case bisC:
    mdat = &bisCMeth;
    sdat = &bisCIsOnPlusStrand;
    break;
  case otherC:
    mdat = &otherCMeth;
    sdat = &otherCIsOnPlusStrand;
    break;
  }
  
  for(int i = 0; i < rposvec.size(); ++i){
    // check whether position already in the map
    if(mdat->find(rposvec[i]) == mdat->end()){ // if not found then insert data
      mdat->insert(make_pair(rposvec[i], protectvec[i]));
      sdat->insert(make_pair(rposvec[i], isonplusstrandvec[i]));
    }
    else if(!isSecond){ // if found and data is not for second mate then overwrite.
      (*mdat)[rposvec[i]] = protectvec[i];
      (*sdat)[rposvec[i]] = isonplusstrandvec[i];
    }
  }
  return(1);
};

// convert data to strings
vector<string > fragData::getStringData() {
  string datStr, strandStr;
  for(int rpos = fragStart; rpos <= fragEnd; ++rpos){
    // hierarchy GCH > WCG > bisC>otherC
    if(gchProtect.find(rpos) != gchProtect.end()){
      datStr += gchProtect[rpos] ? "P" : "p";
      strandStr += gchIsOnPlusStrand[rpos] ? "+" : "-";
    } else if(wcgMeth.find(rpos) != wcgMeth.end()){
      datStr += wcgMeth[rpos] ? "C" : "c";
      strandStr += wcgIsOnPlusStrand[rpos] ? "+" : "-";
    } else if(bisCMeth.find(rpos) != bisCMeth.end()){
      datStr += bisCMeth[rpos] ? "B" : "b";
      strandStr += bisCIsOnPlusStrand[rpos] ? "+" : "-";
    } else if(otherCMeth.find(rpos) != otherCMeth.end()){
      datStr += otherCMeth[rpos] ? "O" : "o";
      strandStr += otherCIsOnPlusStrand[rpos] ? "+" : "-";
    } else{
      datStr += ".";
      strandStr += ".";
    }
  }
  
  vector<string > outd{datStr,strandStr};
  return(outd);
  
};


// print
void fragData::print(){
  vector<string > sdat = getStringData();
  Rcpp::Rcout<<"file prefix="<<prefix<<"; fragStart="<<fragStart<<"; fragEnd="<<fragEnd<<endl;
  Rcpp::Rcout<<"  data:"<<sdat[0]<<endl;
  Rcpp::Rcout<<"strand:"<<sdat[1]<<endl;
};

// 3. function to add cooccurrence counts into provided table
void fragData::addCountsToCoocTable(coocCntTable& cTable){
  map<int, bool >::iterator mend = gchProtect.end();
  for(map<int, bool >::iterator it1 = gchProtect.begin(); it1 != mend; ++it1){
    for(map<int, bool >::iterator it2 = it1; it2 != mend; ++it2){
      if(it1->second==0 && it2->second==0)
        cTable.addCount(abs(it1->first - it2->first),0,1);
      else if(it1->second==0 && it2->second==1)
        cTable.addCount(abs(it1->first - it2->first),1,1);
      else if(it1->second==1 && it2->second==0)
        cTable.addCount(abs(it1->first - it2->first),2,1);
      else if(it1->second==1 && it2->second==1)
        cTable.addCount(abs(it1->first - it2->first),3,1);
    }
  }
  
};


// function to add counts to protection statistics table
void fragData::addCountsToProtectStatsTable(protectStatsTable& pStatsTable){
  
  if(!gchProtect.empty()){
    // count protected gchs
    int protSum = 0;
    for(map<int, bool >::iterator it = gchProtect.begin(); it != gchProtect.end(); ++it){
      protSum += it->second;
    }
    
    // add to pStatsTable
    pStatsTable.addCounts(gchProtect.size(),protSum,1);
  }
  
};



// 2. function to clip data. returns 1 if fragment was clipped or 0 otherwise
int fragData::clipProtectData(const int& clip_until_nbg, // clips from both ends until this number of 0 is met
                    const double& max_protect_frac // keeps fragment if percentage of 1s in a fragment is above this number
){
  if(gchProtect.size() == 0)
    return(0);
  
  // count total protect fraction
  int protSum = 0;
  for(map<int, bool >::iterator it = gchProtect.begin(); it != gchProtect.end(); ++it){
    protSum += it->second;
  }
  
  double protFrac = ((double )protSum)/((double)gchProtect.size());
  if(protFrac >= max_protect_frac) // we keep (nearly) fully protected frags
    return(0);
  
  
  // 1. go forward and remove positions with 1 upto nbg is clip_until_nbg
  
  unordered_set<int > keys_to_remove;
  int nbg = 0;
  for(int pos = fragStart; pos <= fragEnd; ++pos){
    if(gchProtect.find(pos) != gchProtect.end()){
      if(gchProtect[pos] == 0){
        nbg++;
        if(nbg >= clip_until_nbg)
          break;
      } 
      keys_to_remove.insert(pos);
      
    }
  }
  
  // 2. same for backward direction
  nbg = 0;
  for(int pos = fragEnd; pos >= fragStart; --pos){
    if(gchProtect.find(pos) != gchProtect.end()){
      if(gchProtect[pos] == 0){
        nbg++;
        if(nbg >=clip_until_nbg)
          break;
      }
      keys_to_remove.insert(pos);
    }
  }
  
  
  // 3. erase keys
  for(unordered_set<int >::iterator it = keys_to_remove.begin(); it != keys_to_remove.end(); ++it){
    gchProtect.erase(*it);
    gchIsOnPlusStrand.erase(*it);
  }
  
  // return number of positions erased
  return(keys_to_remove.size() > 0 ? 1 : 0);
};


double fragData::get_bisC_meth(int min_size){
  if(bisCMeth.size() < min_size)
    return(0);
  
  // count total fraction of methylation C's in bisC context
  int bisC_sum = 0;
  for(map<int, bool >::iterator it = bisCMeth.begin(); it != bisCMeth.end(); ++it){
    bisC_sum += it->second;
  }
  
  double biscmperc = ((double)bisC_sum)/((double)bisCMeth.size());
  return(biscmperc);
};


// get string representation of fragData
string fragData::getStringEncodingForUniqueness(){
  string enc = prefix + to_string(fragStart) + "-" + to_string(fragEnd)+"|";
  vector<map<int , bool > *> mdatpnts{&gchProtect, &wcgMeth, &bisCMeth, &otherCMeth};
  vector<map<int , bool > *> sdatpnts{&gchIsOnPlusStrand, &wcgIsOnPlusStrand, &bisCIsOnPlusStrand, &otherCIsOnPlusStrand};
  
  for(int i = 0; i < mdatpnts.size(); ++i){
    for(map<int, bool >::iterator it = mdatpnts[i]->begin(); it != mdatpnts[i]->end(); ++it){
      enc += to_string(it->first) + ":" + to_string(it->second) + ":" + to_string((*sdatpnts[i])[it->first]) + ";";
    }
    enc += "|";
  }
  return(enc);
};


