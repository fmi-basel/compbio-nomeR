#include "functions_export_rcpp.hpp"


Rcpp::List fetch_cooc_ctable_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::IntegerVector& max_spac,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax){
  // convert input parameters
  vector<string > infiles_ = Rcpp::as<vector<string > >(infiles);
  string regionChr_ = Rcpp::as<string >(regionChr);
  int regionStart_ = Rcpp::as<int >(regionStart);
  int regionEnd_ = Rcpp::as<int >(regionEnd);
  int max_spac_ = Rcpp::as<int >(max_spac);
  string seqstring_ = Rcpp::as<string >(seqstring);
  int seqStart_ = Rcpp::as<int >(seqStart);
  int seqEnd_ = Rcpp::as<int >(seqEnd);
  
  bool remove_nonunique_ = Rcpp::as<bool >(remove_nonunique);
  int clip_until_nbg_ = Rcpp::as<int >(clip_until_nbg);
  double max_protect_frac_ = Rcpp::as<double >(max_protect_frac);
  double max_bisC_meth_ = Rcpp::as<double >(max_bisC_meth);
  int min_bisC_size_ = Rcpp::as<int >(min_bisC_size);
  int mapqMin_ = Rcpp::as<int >(mapqMin);
  int mapqMax_ = Rcpp::as<int >(mapqMax);
  
  // initialize container for reference sequence
  refSeqInfo refseqInfo(seqstring_, seqStart_,seqEnd_);
  
  // initialize data structure for keeping object pointers
  obj_pnts pnts;
  // store pointer to reference sequence container
  pnts.refseq_info = &refseqInfo;
  // store min and amp MAPQ
  pnts.mapqMin = mapqMin_;
  pnts.mapqMax = mapqMax_;
  
  // initialize container for cooccurrence counts and counters
  coocCntTable cTable(max_spac_);
  int n_fetched = 0, n_nonUnique = 0, n_bisfailed = 0, n_analysed = 0, n_clipped = 0;
  
  
  // for each bam file fetch data and add counts to cooccurrence table 
  for(int i = 0; i < infiles_.size(); ++i){
    // initialize container for fragments in region
    regionData regData(regionChr_,regionStart_,regionEnd_);
    pnts.reg_data = &regData;
    pnts.prefix = "f" + to_string(i + 1)+":";
    // fetch data from bam file
    fetch_data_from_bam(infiles_[i],
                        regionChr_,
                        regionStart_,
                        regionEnd_,
                        &pnts);
    n_fetched += regData.size();
    
    // get qnames on non-unique fragments
    set<string > nonUnfrags;
    if(remove_nonunique_)
      nonUnfrags = regData.getNonUniqueQnames();
    n_nonUnique += nonUnfrags.size();
    
    // get qnames with failed bisulfite conversion
    set<string > failedBisConvfrags;
    if(max_bisC_meth_ < 1.0)
      failedBisConvfrags = regData.getFailedBisConvQnames(max_bisC_meth_,min_bisC_size_);
    n_bisfailed += failedBisConvfrags.size();
    
    // erase fragments which nonunique and failed bisConv
    std::set<string > remove_qnames(nonUnfrags);
    remove_qnames.insert(failedBisConvfrags.begin(), failedBisConvfrags.end());
    n_analysed += regData.eraseQnames(remove_qnames);
    
    // clip fragments
    if(clip_until_nbg_ > 0)
      n_clipped += regData.clipProtectData(clip_until_nbg_,max_protect_frac_);
    
    // add counts to table
    regData.addCountsToCoocTable(cTable);
  }
  
  // convert count table to R matrix and return
  Rcpp::IntegerMatrix cntbl(max_spac_,4);
  vector<string > rownm;
  for(int row = 0; row < max_spac_; ++row){
    for(int col = 0; col < 4; ++col)
      cntbl(row, col) = cTable.get_element(row, col);
    rownm.push_back(to_string(row+1));
    
  }
  Rcpp::CharacterVector colnm = {"N00","N01","N10","N11"};
  
  Rcpp::colnames(cntbl) = colnm;
  Rcpp::rownames(cntbl) = Rcpp::wrap(rownm);
  // create List for R
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("nFragsFetched") = Rcpp::wrap(n_fetched),
                                           Rcpp::Named("nFragsNonUnique") = Rcpp::wrap(n_nonUnique),
                                           Rcpp::Named("nFragsBisFailed") = Rcpp::wrap(n_bisfailed),
                                           Rcpp::Named("nFragsAnalyzed") = Rcpp::wrap(n_analysed),
                                           Rcpp::Named("ClippedUntilNbg") = Rcpp::wrap(clip_until_nbg_),
                                           Rcpp::Named("NotClippedProtectAbove") = Rcpp::wrap(max_protect_frac_),
                                           Rcpp::Named("nFragsClipped") = Rcpp::wrap(n_clipped),
                                           Rcpp::Named("CoocCountTable") = cntbl);
  return(data_out);
}





Rcpp::List fetch_data_matrix_from_bams_cpp(const Rcpp::CharacterVector& whichContext,
                                           const Rcpp::CharacterVector& infiles,
                                           const Rcpp::CharacterVector& regionChr,
                                           const Rcpp::IntegerVector& regionStart,
                                           const Rcpp::IntegerVector& regionEnd,
                                           const Rcpp::CharacterVector& seqstring,
                                           const Rcpp::IntegerVector& seqStart,
                                           const Rcpp::IntegerVector& seqEnd,
                                           const Rcpp::LogicalVector& remove_nonunique,
                                           const Rcpp::IntegerVector& clip_until_nbg,
                                           const Rcpp::NumericVector& max_protect_frac,
                                           const Rcpp::NumericVector& max_bisC_meth,
                                           const Rcpp::IntegerVector& min_bisC_size,
                                           const Rcpp::IntegerVector& mapqMin,
                                           const Rcpp::IntegerVector& mapqMax){
  
  // convert input parameters
  string whichContext_ = Rcpp::as<string>(whichContext);
  fragMapType whichMap;
  if(whichContext_ == "GCH")
    whichMap = GCH;
  else if(whichContext_ == "WCG")
    whichMap = WCG;
  else if(whichContext_ == "bisC")
    whichMap = bisC;
  else if(whichContext_ == "otherC")
    whichMap = otherC;
  else if(whichContext_ == "allC")
    whichMap = allC;
  else
    Rcpp::stop("whichContext must be GCH, WCG, bisC, otherC or allC");
  
  vector<string > infiles_ = Rcpp::as<vector<string > >(infiles);
  string regionChr_ = Rcpp::as<string >(regionChr);
  int regionStart_ = Rcpp::as<int >(regionStart);
  int regionEnd_ = Rcpp::as<int >(regionEnd);
  string seqstring_ = Rcpp::as<string >(seqstring);
  int seqStart_ = Rcpp::as<int >(seqStart);
  int seqEnd_ = Rcpp::as<int >(seqEnd);
  
  bool remove_nonunique_ = Rcpp::as<bool >(remove_nonunique);
  int clip_until_nbg_ = Rcpp::as<int >(clip_until_nbg);
  double max_protect_frac_ = Rcpp::as<double >(max_protect_frac);
  double max_bisC_meth_ = Rcpp::as<double >(max_bisC_meth);
  int min_bisC_size_ = Rcpp::as<int >(min_bisC_size);
  int mapqMin_ = Rcpp::as<int >(mapqMin);
  int mapqMax_ = Rcpp::as<int >(mapqMax);
  
  // initialize container for reference sequence
  refSeqInfo refseqInfo(seqstring_, seqStart_,seqEnd_);
  
  // initialize data structure for keeping object pointers
  obj_pnts pnts;
  // store pointer to reference sequence container
  pnts.refseq_info = &refseqInfo;
  // store min and amp MAPQ
  pnts.mapqMin = mapqMin_;
  pnts.mapqMax = mapqMax_;
  
  // create container for fragments and counters
  regionData regData(regionChr_,regionStart_,regionEnd_);
  pnts.reg_data = &regData;
  int n_fetched = 0, n_nonUnique = 0, n_bisfailed = 0, n_analysed = 0, n_clipped = 0;
  
  // fetch data from all bam files
  for(int i = 0; i < infiles_.size(); ++i){
    pnts.prefix = "f" + to_string(i + 1)+":";
    // fetch data from bam file
    fetch_data_from_bam(infiles_[i],
                        regionChr_,
                        regionStart_,
                        regionEnd_,
                        &pnts);
  }
  
  n_fetched = regData.size();
  
  // get qnames on non-unique fragments
  set<string > nonUnfrags;
  if(remove_nonunique_)
    nonUnfrags = regData.getNonUniqueQnames();
  n_nonUnique = nonUnfrags.size();
  
  // get qnames with failed bisulfite conversion
  set<string > failedBisConvfrags;
  if(max_bisC_meth_ < 1.0)
    failedBisConvfrags = regData.getFailedBisConvQnames(max_bisC_meth_,min_bisC_size_);
  n_bisfailed = failedBisConvfrags.size();
  
  // erase fragments which nonunique and failed bisConv
  std::set<string > remove_qnames(nonUnfrags);
  remove_qnames.insert(failedBisConvfrags.begin(), failedBisConvfrags.end());
  n_analysed = regData.eraseQnames(remove_qnames);
  
  // clip fragments
  if(clip_until_nbg_ > 0)
    n_clipped = regData.clipProtectData(clip_until_nbg_,max_protect_frac_);
  
  
  
  // create List for R
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("nFragsFetched") = Rcpp::wrap(n_fetched),
                                           Rcpp::Named("nFragsNonUnique") = Rcpp::wrap(n_nonUnique),
                                           Rcpp::Named("nFragsBisFailed") = Rcpp::wrap(n_bisfailed),
                                           Rcpp::Named("nFragsAnalyzed") = Rcpp::wrap(n_analysed),
                                           Rcpp::Named("ClippedUntilNbg") = Rcpp::wrap(clip_until_nbg_),
                                           Rcpp::Named("NotClippedProtectAbove") = Rcpp::wrap(max_protect_frac_),
                                           Rcpp::Named("nFragsClipped") = Rcpp::wrap(n_clipped),
                                           Rcpp::Named("DataMatrix") = regData.getDataMatrix(whichMap));
  return(data_out);
}



Rcpp::List fetch_dupl_stats_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                          const Rcpp::CharacterVector& regionChr,
                                          const Rcpp::IntegerVector& regionStart,
                                          const Rcpp::IntegerVector& regionEnd,
                                          const Rcpp::CharacterVector& seqstring,
                                          const Rcpp::IntegerVector& seqStart,
                                          const Rcpp::IntegerVector& seqEnd,
                                          const Rcpp::IntegerVector& mapqMin,
                                          const Rcpp::IntegerVector& mapqMax){
  // convert input parameters
  
  vector<string > infiles_ = Rcpp::as<vector<string > >(infiles);
  string regionChr_ = Rcpp::as<string >(regionChr);
  int regionStart_ = Rcpp::as<int >(regionStart);
  int regionEnd_ = Rcpp::as<int >(regionEnd);
  string seqstring_ = Rcpp::as<string >(seqstring);
  int seqStart_ = Rcpp::as<int >(seqStart);
  int seqEnd_ = Rcpp::as<int >(seqEnd);
  
  int mapqMin_ = Rcpp::as<int >(mapqMin);
  int mapqMax_ = Rcpp::as<int >(mapqMax);
  
  // initialize container for reference sequence
  refSeqInfo refseqInfo(seqstring_, seqStart_,seqEnd_);
  
  // initialize data structure for keeping object pointers
  obj_pnts pnts;
  // store pointer to reference sequence container
  pnts.refseq_info = &refseqInfo;
  // store min and amp MAPQ
  pnts.mapqMin = mapqMin_;
  pnts.mapqMax = mapqMax_;
  
  // create container for fragments and counters
  regionData regData(regionChr_,regionStart_,regionEnd_);
  pnts.reg_data = &regData;
  int n_fetched = 0, n_nonUnique = 0;
  
  // fetch data from all bam files
  for(int i = 0; i < infiles_.size(); ++i){
    pnts.prefix = "f" + to_string(i + 1)+":";
    // fetch data from bam file
    fetch_data_from_bam(infiles_[i],
                        regionChr_,
                        regionStart_,
                        regionEnd_,
                        &pnts);
  }
  
  n_fetched = regData.size();
  
  // get qnames on non-unique fragments. this is not essential in this function, but for testing is useful
  set<string > nonUnfrags;
  nonUnfrags = regData.getNonUniqueQnames();
  n_nonUnique = nonUnfrags.size();
  
  // get map duplication encoding -> vector of qnames
  unordered_map<string, vector<string > > duplEncToQnamesMap = regData.getDuplQnames();
  
  // get stats of how many each duplication encoding occurs
  map<int, int> dupl_stats; // map dupl_times -> number of fragments duplicated this number of times
  for(unordered_map<string, vector<string > >::iterator it = duplEncToQnamesMap.begin(); it != duplEncToQnamesMap.end(); ++it){
    // get number of frags for current encoding
    int nfrgs = (it->second).size();
    
    // check if nfrgs exists in dupl_stats and add if not
    if(dupl_stats.find(nfrgs) == dupl_stats.end()){ // this number of duplications is not found
      // insert a new pair nfrgs -> 1
      dupl_stats.insert(make_pair(nfrgs,1));
    } else{ // this number of duplications is already in the map
      dupl_stats[nfrgs]++;
    }
  }
  
  // create a matrix for duplication stats: 
  // columns
  // number_of_frags -> n_of_times duplicated
  Rcpp::IntegerMatrix duplStats4R(dupl_stats.size(),2);
  fill(duplStats4R.begin(),duplStats4R.end(),NA_INTEGER);
  Rcpp::CharacterVector colnms = {"n_frags","n_times_duplicated"};
  Rcpp::colnames(duplStats4R) = colnms;
  int rowi=0; // row index in the matrix
  for(map<int, int>::iterator it = dupl_stats.begin(); it != dupl_stats.end(); ++it){
    duplStats4R(rowi,0) = it->first;
    duplStats4R(rowi,1) = it->second;
    ++rowi;
  }
  
  // create Rcpp List for returning unique encodinds and fragnames
  Rcpp::List EncToQnames4R;
  for(unordered_map<string, vector<string > >::iterator it = duplEncToQnamesMap.begin(); it != duplEncToQnamesMap.end(); ++it){
    EncToQnames4R.push_back(Rcpp::wrap(it->second),it->first);
  }
  // create List for R
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("nFragsFetched") = Rcpp::wrap(n_fetched),
                                           Rcpp::Named("nFragsNonUnique") = Rcpp::wrap(n_nonUnique),
                                           Rcpp::Named("duplStats") = Rcpp::wrap(duplStats4R),
                                           Rcpp::Named("duplFragNames") = EncToQnames4R);
  return(data_out);
}




Rcpp::List fetch_protect_stats_from_bams_cpp(const Rcpp::CharacterVector& infiles,
                                             const Rcpp::CharacterVector& regionChr,
                                             const Rcpp::IntegerVector& regionStart,
                                             const Rcpp::IntegerVector& regionEnd,
                                             const Rcpp::CharacterVector& seqstring,
                                             const Rcpp::IntegerVector& seqStart,
                                             const Rcpp::IntegerVector& seqEnd,
                                             const Rcpp::LogicalVector& remove_nonunique,
                                             const Rcpp::IntegerVector& clip_until_nbg,
                                             const Rcpp::NumericVector& max_protect_frac,
                                             const Rcpp::NumericVector& max_bisC_meth,
                                             const Rcpp::IntegerVector& min_bisC_size,
                                             const Rcpp::IntegerVector& mapqMin,
                                             const Rcpp::IntegerVector& mapqMax){
  
  // convert input parameters
  vector<string > infiles_ = Rcpp::as<vector<string > >(infiles);
  string regionChr_ = Rcpp::as<string >(regionChr);
  int regionStart_ = Rcpp::as<int >(regionStart);
  int regionEnd_ = Rcpp::as<int >(regionEnd);
  string seqstring_ = Rcpp::as<string >(seqstring);
  int seqStart_ = Rcpp::as<int >(seqStart);
  int seqEnd_ = Rcpp::as<int >(seqEnd);
  
  bool remove_nonunique_ = Rcpp::as<bool >(remove_nonunique);
  int clip_until_nbg_ = Rcpp::as<int >(clip_until_nbg);
  double max_protect_frac_ = Rcpp::as<double >(max_protect_frac);
  double max_bisC_meth_ = Rcpp::as<double >(max_bisC_meth);
  int min_bisC_size_ = Rcpp::as<int >(min_bisC_size);
  int mapqMin_ = Rcpp::as<int >(mapqMin);
  int mapqMax_ = Rcpp::as<int >(mapqMax);
  
  
  // initialize container for reference sequence
  refSeqInfo refseqInfo(seqstring_, seqStart_,seqEnd_);
  
  // initialize data structure for keeping object pointers
  obj_pnts pnts;
  // store pointer to reference sequence container
  pnts.refseq_info = &refseqInfo;
  // store min and max MAPQ
  pnts.mapqMin = mapqMin_;
  pnts.mapqMax = mapqMax_;
  
  // initialize container for collecting protect statistics and counters
  protectStatsTable protStats;
  int n_fetched = 0, n_nonUnique = 0, n_bisfailed = 0, n_analysed = 0, n_clipped = 0;
  
  
  // 	// for each bam file fetch data and add counts to protect stats table 
  for(int i = 0; i < infiles_.size(); ++i){
    // initialize container for fragments in region
    regionData regData(regionChr_,regionStart_,regionEnd_);
    pnts.reg_data = &regData;
    pnts.prefix = "f" + to_string(i + 1)+":";
    // fetch data from bam file
    fetch_data_from_bam(infiles_[i],
                        regionChr_,
                        regionStart_,
                        regionEnd_,
                        &pnts);
    n_fetched += regData.size();
    
    // get qnames on non-unique fragments
    set<string > nonUnfrags;
    if(remove_nonunique_)
      nonUnfrags = regData.getNonUniqueQnames();
    n_nonUnique += nonUnfrags.size();
    
    // get qnames with failed bisulfite conversion
    set<string > failedBisConvfrags;
    if(max_bisC_meth_ < 1.0)
      failedBisConvfrags = regData.getFailedBisConvQnames(max_bisC_meth_,min_bisC_size_);
    n_bisfailed += failedBisConvfrags.size();
    
    // erase fragments which nonunique and failed bisConv
    std::set<string > remove_qnames(nonUnfrags);
    remove_qnames.insert(failedBisConvfrags.begin(), failedBisConvfrags.end());
    n_analysed += regData.eraseQnames(remove_qnames);
    
    // clip fragments
    if(clip_until_nbg_ > 0)
      n_clipped += regData.clipProtectData(clip_until_nbg_,max_protect_frac_);
    
    // add counts to protect stat table
    regData.addCountsToProtectStatsTable(protStats);
  }
  
  // convert stat table to R matrix
  vector<vector<int >> stTbl = protStats.getStatTable();
  
  Rcpp::IntegerMatrix pStatMatr(stTbl.size(),3);
  for(int i=0;i<stTbl.size(); ++i){
    pStatMatr(i,0) = stTbl[i][0];
    pStatMatr(i,1) = stTbl[i][1];
    pStatMatr(i,2) = stTbl[i][2];
  }
  
  Rcpp::CharacterVector colnm = {"TotalGCH","ProtectedGCH","Nfrags"};
  
  Rcpp::colnames(pStatMatr) = colnm;
  
  // create List for R
  Rcpp::List data_out = Rcpp::List::create(Rcpp::Named("nFragsFetched") = Rcpp::wrap(n_fetched),
                                           Rcpp::Named("nFragsNonUnique") = Rcpp::wrap(n_nonUnique),
                                           Rcpp::Named("nFragsBisFailed") = Rcpp::wrap(n_bisfailed),
                                           Rcpp::Named("nFragsAnalyzed") = Rcpp::wrap(n_analysed),
                                           Rcpp::Named("ClippedUntilNbg") = Rcpp::wrap(clip_until_nbg_),
                                           Rcpp::Named("NotClippedProtectAbove") = Rcpp::wrap(max_protect_frac_),
                                           Rcpp::Named("nFragsClipped") = Rcpp::wrap(n_clipped),
                                           Rcpp::Named("ProtectStats") = pStatMatr);
  return(data_out);
}




Rcpp::List run_cpp_nomeR(const Rcpp::List& data,
                         const Rcpp::CharacterVector& fragnames,
                         const Rcpp::List& binding_models,
                         const Rcpp::NumericVector& bgprotectprob,
                         const Rcpp::NumericVector& bgprior,
                         const Rcpp::LogicalVector& report_prediction_in_flanks,
                         const Rcpp::NumericVector& Ncpu,
                         const Rcpp::LogicalVector& verbose
) {
  
  //set verbose
  extern bool _VERBOSE_;
  _VERBOSE_ = Rcpp::as<bool >(verbose);
  
  // set report_prediction_in_flanks
  bool report_prediction_in_flanks_ = Rcpp::as<bool >(report_prediction_in_flanks);
  
  int Ncpu_ = Rcpp::as<int >(Ncpu);
#ifndef _OPENMP
  Rcpp::Rcout<<"nomeR was compiled without OpenMP. ncpu does not have effect.\n";
#endif
  
  
  
  if(_VERBOSE_){
    Rcpp::Rcout<<"Creating Predict object..."<<endl;
  }
  Predict predict(data,
                  fragnames,
                  binding_models,
                  bgprotectprob,
                  bgprior);
  
  // run prediction
  if(_VERBOSE_){
    Rcpp::Rcout<<"Calculating posterior binding probabilities..."<<endl;
  }
  
  Rcpp::List output_data;
  if(predict.Run(Ncpu_)){
    //predict.Run(Ncpu_);
    
    if(_VERBOSE_){
      Rcpp::Rcout<<"Running predict.getStartProbDF()..."<<endl;
    }
    Rcpp::List startProbdf = predict.getStartProbDF(report_prediction_in_flanks_);
    if(_VERBOSE_){
      Rcpp::Rcout<<"Running predict.getCoverProbDF()..."<<endl;
    }
    Rcpp::List coverProb = predict.getCoverProbDF();
    if(_VERBOSE_){
      Rcpp::Rcout<<"Running predict.getGenomeSummaryDF()..."<<endl;
    }
    Rcpp::List genSummary = predict.getGenomeSummaryDF();
    output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = startProbdf,
                                      Rcpp::Named("COVER_PROB") = coverProb,
                                      Rcpp::Named("SUMMARY") = genSummary);
  } else {
    output_data = Rcpp::List::create( Rcpp::Named("START_PROB") = R_NilValue,
                                      Rcpp::Named("COVER_PROB") = R_NilValue,
                                      Rcpp::Named("SUMMARY") = R_NilValue);
    
  }
  
  _VERBOSE_ = 0;
  return output_data;
  
}


Rcpp::List count_spacing_freq_cpp(const Rcpp::List& data,
                            const Rcpp::CharacterVector& fragnames,
                            const Rcpp::IntegerVector& maxspacing,
                            const Rcpp::IntegerVector& maxwmlen){
  
  int maxwmlen_ = Rcpp::as<int >(maxwmlen);
  int maxspacing_ = Rcpp::as<int >(maxspacing);
  
  NOMeSeqData nome_data(data,
                        fragnames,
                        maxwmlen_);
  //nome_data.create(data,maxwmlen_);
  vector<vector<int> > freq_mat = nome_data.count_freq_for_spacings(maxspacing_);
  vector<int > spacings;
  vector<int > freq00;
  vector<int > freq01;
  vector<int > freq02;
  vector<int > freq10;
  vector<int > freq11;
  vector<int > freq12;
  
  vector<int > freq20;
  vector<int > freq21;
  vector<int > freq22;
  for(int i=0; i<freq_mat.size(); ++i){
    spacings.push_back(freq_mat[i][0]);
    freq00.push_back(freq_mat[i][1]);
    freq01.push_back(freq_mat[i][2]);
    freq02.push_back(freq_mat[i][3]);
    freq10.push_back(freq_mat[i][4]);
    freq11.push_back(freq_mat[i][5]);
    freq12.push_back(freq_mat[i][6]);
    freq20.push_back(freq_mat[i][7]);
    freq21.push_back(freq_mat[i][8]);
    freq22.push_back(freq_mat[i][9]);
    
  }
  
  Rcpp::List export_list;
  
  export_list.push_back(Rcpp::wrap(spacings),"S");
  export_list.push_back(Rcpp::wrap(freq00),"N00");
  export_list.push_back(Rcpp::wrap(freq01),"N01");
  export_list.push_back(Rcpp::wrap(freq10),"N10");
  export_list.push_back(Rcpp::wrap(freq11),"N11");
  export_list.push_back(Rcpp::wrap(freq02),"N0na");
  export_list.push_back(Rcpp::wrap(freq12),"N1na");
  export_list.push_back(Rcpp::wrap(freq20),"Nna0");
  export_list.push_back(Rcpp::wrap(freq21),"Nna1");
  export_list.push_back(Rcpp::wrap(freq22),"Nnana");
  
  
  nome_data.clear();
  return export_list;
  
}


Rcpp::List calculate_theor_joint_prob_cpp(const Rcpp::NumericVector& ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1. Make sure that R function passes correct vector with priors
                                    const Rcpp::NumericVector& bg_protect_prob,
                                    const Rcpp::NumericVector& footprint_protect_prob,
                                    const Rcpp::IntegerVector& max_spacing){
  
  vector<double > ftp_cover_priors_ = Rcpp::as<vector<double > >(ftp_cover_priors);
  double bg_protect_prob_ = Rcpp::as<double >(bg_protect_prob);
  double footprint_protect_prob_ = Rcpp::as<double >(footprint_protect_prob);
  int max_spacing_ = Rcpp::as<int >(max_spacing);
  
  
  DNAbind_obj_vector ftp_array;
  vector<vector<double > > p_joint = ftp_array.calc_theor_joint_prob(ftp_cover_priors_,
                                                                     bg_protect_prob_,
                                                                     footprint_protect_prob_,
                                                                     max_spacing_);
  
  // create Rcpp object
  
  vector<int > spacings;
  vector<double > p00;
  vector<double > p01;
  vector<double > p10;
  vector<double > p11;
  
  
  for(int i = 0; i < p_joint.size(); ++i){
    //Rcpp::Rcout<<"S="<<p_joint[i][0]<<"; "<<p_joint[i][1]<<"; "<<p_joint[i][2]<<"; "<<p_joint[i][3]<<"; "<<p_joint[i][4]<<endl;
    spacings.push_back(p_joint[i][0]);
    p00.push_back(p_joint[i][1]);
    p01.push_back(p_joint[i][2]);
    p10.push_back(p_joint[i][3]);
    p11.push_back(p_joint[i][4]);
  }
  
  Rcpp::List export_list;
  export_list.push_back(Rcpp::wrap(spacings),"S");
  export_list.push_back(Rcpp::wrap(p00),"P00");
  export_list.push_back(Rcpp::wrap(p01),"P01");
  export_list.push_back(Rcpp::wrap(p10),"P10");
  export_list.push_back(Rcpp::wrap(p11),"P11");
  
  
  ftp_array.clear();
  return export_list;
}




