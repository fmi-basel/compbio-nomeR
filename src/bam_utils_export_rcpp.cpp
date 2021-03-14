#include "bam_utils_export_rcpp.hpp"


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


