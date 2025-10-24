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
                       const Rcpp::IntegerVector& protectVec,
                       int maxWMlen){
	create(fragIDs,
        fragPos,
        protectVec,
        maxWMlen);

}

bool SMFdataset::create(const Rcpp::IntegerVector& fragIDs,
                        const Rcpp::IntegerVector& fragPos,
                        const Rcpp::IntegerVector& protectVec,
                        int maxWMlen){
	_nmolecs = 0;
	_totallength = 0;
	if(fragIDs.size() != fragPos.size() || fragIDs.size() != protectVec.size())
		Rcpp::stop("SMFdataset::create: Inconsistent lengths of input vectors fragIDs, fragPos and protecVec.\n");

	int currFragID = -1;
	vector<uint32_t > currfragPosVec; // input positions fragPosVec must be 0-based
	vector<uint8_t > currprotectVec;
	for(int i = 0; i < fragIDs.size(); ++i){

		uint32_t fragid = fragIDs[i];
		uint32_t fragpos = fragPos[i];
		uint8_t protectval = protectVec[i];

		if(fragid != currFragID){
			if(currFragID == -1){ // beginning of the loop
				currFragID = fragid;
			} else{ // beginning of data for a new fragment
				// add old fragment
				Add(currFragID,currfragPosVec, currprotectVec, maxWMlen);
				// clear vectors
				currfragPosVec.clear();
				currprotectVec.clear();
				// assign new curFragID
				currFragID = fragid;
			}

		}

		// add fragpos and protectval
		currfragPosVec.push_back(fragpos);
		currprotectVec.push_back(protectval);
	}
	// add the last fragment
	Add(currFragID,currfragPosVec, currprotectVec, maxWMlen);

	return 1;
}

void SMFdataset::Add(const uint32_t fragID,
                     const vector<uint32_t>& fragPosVec, // input positions fragPosVec must be 0-based
                     const vector<uint8_t>& protectVec,
                     int maxWMlen){
	// create a new object of class fragProtectData
	fragProtectData newFrag(fragID,fragPosVec, protectVec, maxWMlen);
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

fragProtectData & SMFdataset::operator[] (int index)
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


vector<vector<int> > SMFdataset::count_freq_for_spacings(int maxSpacing) const
{
	// here the spacing 0 means that positions are adjacent and gap between them is 0
	// however output S will be starting from 1
	// create output vector of vectors
	vector<vector<int > > freq_mat;
	// for each spacing 0:maxSpacing
	for(int s=0; s < maxSpacing; ++s){
		vector<int > freq_vec(5,0); // columns are S; 0,0; 0,1; 1,0; 1,1;
		freq_vec[0] = s + 1;
		// for each sequence
		for(int f = 0; f < _nmolecs; ++f){
			fragProtectData fragData = _data[f];
			uint32_t firstDatPos = fragData._firstDatpos;
			uint32_t lastDatPos = fragData._lastDatpos;
			// go from first position to the last - s + 1
			for(int pos = firstDatPos; pos <= lastDatPos - s; pos++){
				int letter_pos = fragData[pos];
				int letter_spac = fragData[pos + s];

				if(letter_pos == 0 && letter_spac == 0){
					freq_vec[1]++;
				} else if(letter_pos == 0 && letter_spac == 1){
					freq_vec[2]++;
				} else if(letter_pos == 1 && letter_spac == 0){
					freq_vec[4]++;
				} else if(letter_pos == 1 && letter_spac == 1){
					freq_vec[5]++;
				}

			}
		}

		freq_mat.push_back(freq_vec);
	}
	return freq_mat;
}

