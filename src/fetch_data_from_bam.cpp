#include "fetch_data_from_bam.hpp"

// fetch data from bam file
bool fetch_data_from_bam(const string& bamfile,
                         const string& regionChr,
                         const int& regionStart,
                         const int& regionEnd,
                         obj_pnts* pnts){
	
	
	bam1_t *hit = bam_init1();
	samfile_t *fin;
	bam_index_t *idx;
	
	
	fin = _bam_tryopen(bamfile.c_str(), "rb", NULL);
	idx = bam_index_load(bamfile.c_str()); // load BAM index
	if (idx == 0)
		Rcpp::stop("BAM index for '" + bamfile + "' unavailable\n");
	
	// get target id
	int tid = 0;
	while(regionChr.compare(fin->header->target_name[tid]) != 0 && tid+1<fin->header->n_targets)
		tid++;
	
	if(regionChr.compare(fin->header->target_name[tid]) != 0)
		Rcpp::stop("could not find target '" + regionChr + "' in bam header of " + bamfile + ".\n");
	
	// call collectRegionData on store alignments in region
	bam_fetch(fin->x.bam, idx, tid, regionStart, regionEnd, pnts, &collectRegionData);
	
	// clean bam file objects
	bam_index_destroy(idx);
	samclose(fin);
	
	return(1);
}


static int collectRegionData(const bam1_t *hit, void *data) { // bam_fetch callback function that parses XM tag in bismark bam and adds cooccurrence counts
	
	static obj_pnts *pdata = NULL;
	
	pdata = (obj_pnts*) data;
	
	
	// check if alignment is valid
	if(!isValidAln(hit,pdata))
		return(0);
	// for testing only
	//print_bs_alignment(hit,pdata);
	
	// parse alignment and create vectors for regData
	string qname = string(bam1_qname(hit)); // frag name
	uint32_t *cigar = bam_get_cigar(hit); // cigar array
	uint8_t *hitseq = bam1_seq(hit); // query sequence
	int ref_start = hit -> core.pos; // leftmost position in reference
	int ref_end = bam_calend(&(hit->core),cigar); // rightmost position in reference
	bool is_read_rev = bam_is_rev(hit); // is read on "+" strand
	bool is_mate_rev = bam_is_mrev(hit);
	bool is_second = bam_is_second(hit); // is it a second mate
	
	// choose strand to get contexts for
	seekMismType msmType;
	if(is_second){ // parsing second mate in PE data, strand determined by first mate
		msmType = is_mate_rev ? GtoA : CtoT;
	} else{ // parsing first mate or SE data, strand determined by read itself
		msmType = is_read_rev ? GtoA : CtoT;
	}
	// Rcpp::Rcout<<"msmType="<<msmType<<endl;
	// 
	// // for testing only
	// fragData frg;
	
	// go across contexts and add data
	vector<fragMapType > context_vec = {GCH, WCG, bisC, otherC};
	for(int i = 0; i < context_vec.size(); ++i){
		// get positions for the context in reference
		vector<int> cntx_rpos_vec = pdata->refseq_info->getPositionsBetween(ref_start,ref_end,context_vec[i], msmType == CtoT);
		
		vector<int > rpos_vec;
		vector<bool > dat_vec;
		vector<bool > isplusstrand_vec;
		
		for(int cind = 0; cind < cntx_rpos_vec.size(); ++cind){
			
			// get query position
			int qpos = _fromRefToQueryPos(ref_start,
                                 cntx_rpos_vec[cind],
                                              hit->core.n_cigar,
                                              cigar);
			// get nucleotide in query
			string qseq = nt16_table[bam1_seqi(hitseq, qpos)];
			
			// check if C->T mismatch if msmType == CtoT or GtoA otherwise. query on - in BAM are reverse complement
			switch(msmType) {
				case CtoT:
				{
					if(qseq == "C"){ // methylated or unprotected
					// add position on reference
					rpos_vec.push_back(cntx_rpos_vec[cind]);
					// add data
					dat_vec.push_back(context_vec[i] == GCH ? 0 : 1);
					// add strand
					isplusstrand_vec.push_back(!is_read_rev);
					
				} else if(qseq == "T"){ // unmethylated or protected
					// add position on reference
					rpos_vec.push_back(cntx_rpos_vec[cind]);
					// add data
					dat_vec.push_back(context_vec[i] == GCH ? 1 : 0);
					// add strand
					isplusstrand_vec.push_back(!is_read_rev);
				}
					break;
				}
				case GtoA:
				{
					if(qseq == "G"){ // methylated or unprotected
					// add position on reference
					rpos_vec.push_back(cntx_rpos_vec[cind]);
					// add data
					dat_vec.push_back(context_vec[i] == GCH ? 0 : 1);
					// add strand
					isplusstrand_vec.push_back(!is_read_rev);
					
				} else if(qseq == "A"){ // unmethylated or protected
					// add position on reference
					rpos_vec.push_back(cntx_rpos_vec[cind]);
					// add data
					dat_vec.push_back(context_vec[i] == GCH ? 1 : 0);
					// add strand
					isplusstrand_vec.push_back(!is_read_rev);
				}
					
					break;
				}
			}
			
			
		}
		
		
			
		
		// add data for the current context into regionData
		pdata->reg_data->addFragContextData(pdata->prefix,
                                      pdata->prefix + qname,
                                      ref_start,
                                      ref_end,
                                      rpos_vec,
                                      dat_vec,
                                      isplusstrand_vec,
                                      is_second,
                                      context_vec[i]);
		
		// for testing only
// 		frg.addData(ref_start,
//               ref_end,
//               rpos_vec,
//               dat_vec,
//               isplusstrand_vec,
//               is_second,
//               context_vec[i]);
		
	}
	
	
	
	

	
	return(0);
	
}


// checks validity of an alignment
bool isValidAln(const bam1_t *hit,const obj_pnts* pdata){
	
	return(!bam_is_not_primary(hit) && 
         !bam_is_qc_fail(hit) && 
         !bam_is_pcr_dupl(hit) &&
         !bam_is_suppl(hit) &&
         ((int )(hit->core).qual >= pdata->mapqMin) &&
         ((int )(hit->core).qual <= pdata->mapqMax)
        );
}


// print alignment for debugging
void print_bs_alignment(const bam1_t *hit,const obj_pnts* pdata){
	
	// print data for one alignment
	Rcpp::Rcout<<"====== "<<bam1_qname(hit)<<" ======"<<endl;
	Rcpp::Rcout<<"alignment is first: "<<bam_is_first(hit)<<"; alignment is reverse: "<<bam_is_rev(hit)<<endl;
	
	bool is_read_rev = bam_is_rev(hit); // is read on "+" strand
	bool is_mate_rev = bam_is_mrev(hit);
	bool is_second = bam_is_second(hit); // is it a second mate
	
	// get cigar array
	uint32_t *cigar = bam_get_cigar(hit);
	
	// get qseq, rseq and xm tag and print out
  // old: int qlen = bam_cigar2qlen(&(hit -> core), cigar);
  int qlen = bam_cigar2qlen(hit->core.n_cigar, cigar);
	uint8_t *hitseq = bam1_seq(hit);
	
	
	// ## for Bismark meth and unmeth alphabet find positions in XM tag
	// ## Bismark codes
	// # z - C in CpG context - unmethylated
	// # Z - C in CpG context - methylated
	// # x - C in CHG context - unmethylated
	// # X - C in CHG context - methylated
	// # h - C in CHH context - unmethylated
	// # H - C in CHH context - methylated
	// # u - C in Unknown context (CN or CHN) - unmethylated
	// # U - C in Unknown context (CN or CHN) - methylated
	// # . - not a C or irrelevant position
	
	
	// get XM tag if present
	uint8_t *pxmtag = bam_aux_get(hit,"XM");
	string xm_tag;
	if(pxmtag != NULL){
		int type = *pxmtag++;
		if(type == 'Z')
			xm_tag = string((char*)pxmtag);
	}
	
	
	// choose strand to get contexts for
	seekMismType msmType;
	if(is_second){ // parsing second mate in PE data, strand determined by first mate
		msmType = is_mate_rev ? GtoA : CtoT;
	} else{ // parsing first mate or SE data, strand determined by read itself
		msmType = is_read_rev ? GtoA : CtoT;
	}
	
	string queryseq, refseq, matchsymb,cntxstr;
	for(int qpos = 0; qpos < qlen; ++qpos){
		int rpos = _fromQueryToRefPos(hit -> core.pos,qpos,hit -> core.n_cigar,cigar);
		string rseq;
		string qseq = nt16_table[bam1_seqi(hitseq, qpos)];
		
		if(rpos != -1){
			rseq = pdata->refseq_info->getRefSubseq(rpos,1);
			matchsymb += "|";
		} else {
			rseq = "-";
			matchsymb += " ";
		}
		
		queryseq += qseq;
		refseq += rseq;
		
		// get context
		short int mcntx = pdata->refseq_info->getContextForRefPos(rpos);
		
		switch (msmType) {
		case CtoT:
		{
			if(mcntx == 0 || (qseq != "C" && qseq != "T")){
				cntxstr += ".";
			} else{
				if(rseq != "C"){
					Rcpp::Rcout<<"print_bs_alignment: WARNING: Expect C in reference sequence but got "<<rseq<<endl;
					cntxstr += "X";
				} else if(mcntx == 1){
					if(qseq == "T")
						cntxstr += "P";
					else if(qseq == "C")
						cntxstr += "p";
				} else if(mcntx == 2){
					if(qseq == "T")
						cntxstr += "c";
					else if(qseq == "C")
						cntxstr += "C";
				} else if(mcntx == 3){
					if(qseq == "T")
						cntxstr += "b";
					else if(qseq == "C")
						cntxstr += "B";
				} else if(mcntx == 4){
					if(qseq == "T")
						cntxstr += "o";
					else if(qseq == "C")
						cntxstr += "O";
				}
			
			}
			break;
		}
		case GtoA:
		{
			if(mcntx == 0 || (qseq != "G" && qseq != "A")){
			cntxstr += ".";
		} else{
			if(rseq != "G"){
				Rcpp::Rcout<<"print_bs_alignment: WARNING: Expect G in reference sequence but got "<<rseq<<endl;
				cntxstr += "X";
			} else if(mcntx == -1){
				if(qseq == "A")
					cntxstr += "P";
				else if(qseq == "G")
					cntxstr += "p";
			} else if(mcntx == -2){
				if(qseq == "A")
					cntxstr += "c";
				else if(qseq == "G")
					cntxstr += "C";
			} else if(mcntx == -3){
				if(qseq == "A")
					cntxstr += "b";
				else if(qseq == "G")
					cntxstr += "B";
			} else if(mcntx == -4){
				if(qseq == "A")
					cntxstr += "o";
				else if(qseq == "G")
					cntxstr += "O";
			}
			
		}
			
			break;
		}
		}
		
		
		
	}
	
	// print reference sequence
	Rcpp::Rcout<<"reference:  "<<refseq<<endl;
	// print match symbols
	Rcpp::Rcout<<"            "<<matchsymb<<endl;
	// print query seq
	Rcpp::Rcout<<"    query:  "<<queryseq<<endl;
	// print nome context
	Rcpp::Rcout<<"  context:  "<<cntxstr<<endl;
	
	// print XM tag if present
	if(xm_tag != "")
		Rcpp::Rcout<<"   XM tag:  "<<xm_tag<<endl;
	// print end of alignment 
	Rcpp::Rcout<<"======================================================================================="<<endl;
}

