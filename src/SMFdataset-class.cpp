#include "SMFdataset-class.hpp"

SMFdataset::SMFdataset()
{
	_nmolecs = 0;
	_totallength = 0;
}
SMFdataset::~SMFdataset()
{

}

SMFdataset::SMFdataset(const Rcpp::IntegerVector& fragIDs,
                       const Rcpp::IntegerVector& fragPos,
                       const Rcpp::NumericVector& modProbVec,
                       int maxWMlen){
	create(fragIDs,
        fragPos,
        modProbVec,
        maxWMlen);

}

bool SMFdataset::create(const Rcpp::IntegerVector& fragIDs,
                        const Rcpp::IntegerVector& fragPos,
                        const Rcpp::NumericVector& modProbVec,
                        int maxWMlen){
	_nmolecs = 0;
	_totallength = 0;
	if(fragIDs.size() != fragPos.size() || fragIDs.size() != modProbVec.size())
		Rcpp::stop("SMFdataset::create: Inconsistent lengths of input vectors fragIDs, fragPos and modProbVec.\n");

	int currFragID = -1;
	vector<uint32_t > currfragPosVec; // input positions fragPosVec must be 1-based
	vector<double > currModProbVec;
	for(int i = 0; i < fragIDs.size(); ++i){

		uint32_t fragid = fragIDs[i];
		uint32_t fragpos = fragPos[i];
		double modprobval = modProbVec[i];

		if(fragid != currFragID){
			if(currFragID == -1){ // beginning of the loop
				currFragID = fragid;
			} else{ // beginning of data for a new fragment
				// add old fragment
				Add(currFragID, currfragPosVec, currModProbVec, maxWMlen);
				// clear vectors
				currfragPosVec.clear();
				currModProbVec.clear();
				// assign new curFragID
				currFragID = fragid;
			}

		}

		// add fragpos and modprobval
		currfragPosVec.push_back(fragpos);
		currModProbVec.push_back(modprobval);
	}
	// add the last fragment
	Add(currFragID, currfragPosVec, currModProbVec, maxWMlen);

	return 1;
}

void SMFdataset::Add(const uint32_t fragID,
                     const vector<uint32_t>& fragPosVec, // input positions fragPosVec must be 1-based
                     const vector<double>& modProbVec,
                     int maxWMlen){
	// create a new object of class fragProtectData
	fragProtectData newFrag(fragID, fragPosVec, modProbVec, maxWMlen);
	_totallength += newFrag.Size();
	_data.push_back(newFrag);
	_nmolecs++;
}

void SMFdataset::Add(fragProtectData & frag){
	fragProtectData newFrag(frag);
	_totallength += newFrag.Size();
	_data.push_back(newFrag);
	_nmolecs++;
}

void SMFdataset::clear(){
	_data.clear();
	_nmolecs = 0 ;
	_totallength = 0;
}

const fragProtectData & SMFdataset::operator[] (int index) const
{
	return _data[index];
}

int SMFdataset::Size() const
{
	return _nmolecs;
}

int SMFdataset::TotalLength() const
{
	return _totallength;
}


vector<vector<double > > SMFdataset::count_freq_for_spacings(int maxSpacing,
                                                             int ncpu) const
{
	extern bool _VERBOSE_;
	// spacing 0 means adjacent positions; output S starts from 1
	// columns: E[N(acc,acc)], E[N(acc,prot)], E[N(prot,acc)], E[N(prot,prot)]
	// where acc = accessible (high mod_prob p), prot = protected (low mod_prob, 1-p)
	vector<vector<double > > freqM_glob(maxSpacing, vector<double>(4, 0.0));


#ifdef _OPENMP
	omp_set_nested(true);
	omp_set_num_threads(ncpu);
	if(_VERBOSE_)
		Rcpp::Rcout<<"Running aggregation of co-occurrence statistics with "<<omp_get_max_threads()<<" cpu."<<endl;
#endif

	// parallelize for each molecule
	int seq = 0;
#pragma omp parallel private(seq)
{
	// Each thread gets a private local matrix
	vector<vector<double > > freqM_loc(maxSpacing, vector<double>(4, 0.0));
#pragma omp for schedule(dynamic)
	// for each sequence
	for(seq = 0; seq < _nmolecs; ++seq){
		fragProtectData fragData = _data[seq];
		uint32_t firstDatPos = fragData._firstDatpos;
		uint32_t lastDatPos = fragData._lastDatpos;
		// for each spacing
		for(int s = 0; s < maxSpacing; ++s){
			// go from first position to the last - s + 1
			for(int pos = firstDatPos; (int)(pos + s) <= (int)lastDatPos; ++pos){
				double p_i = fragData[pos];
				double p_j = fragData[pos + s];
				// skip NA-encoded positions (sentinel -1.0)
				if(p_i < 0.0 || p_j < 0.0) continue;
				freqM_loc[s][0] += p_i         * p_j;         // E[N(accessible, accessible)]
				freqM_loc[s][1] += p_i         * (1.0 - p_j); // E[N(accessible, protected)]
				freqM_loc[s][2] += (1.0 - p_i) * p_j;         // E[N(protected, accessible)]
				freqM_loc[s][3] += (1.0 - p_i) * (1.0 - p_j); // E[N(protected, protected)]
			}
		}
	}

	// Reduction: safely combine thread-local matrices into global freqM_glob
#pragma omp critical
{
	for (int s = 0; s < maxSpacing; ++s)
		for (int k = 0; k < 4; ++k)
			freqM_glob[s][k] += freqM_loc[s][k];
}
}

return freqM_glob;
}

