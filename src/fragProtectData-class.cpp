#include "fragProtectData-class.hpp"



fragProtectData::fragProtectData(){

}

fragProtectData::~fragProtectData(){

}

fragProtectData::fragProtectData(const fragProtectData & s){

	_size = s.Size();
	_fragID = s._fragID;
	_firstDatpos = s._firstDatpos;
	_lastDatpos = s._lastDatpos;
	_modProbVec = s._modProbVec;
}

fragProtectData::fragProtectData(const uint32_t fragID,
                                 const vector<uint32_t>& fragPosVec, // fragPosVec are 1 - based
                                 const vector<double>& modProbVec,
                                 int maxWMlen){
	_fragID = fragID;
	// here, we have to redefine positions within the fragment to take into account padding by NAs of size maxWMlen
	// for Posterior-Viterbi decoding we need coverage posteriors in the left flanking region of size maxWMlen
	// to calculate these coverage posteriors we need starting posteriors in the left flank of size 2*maxWMlen
	uint32_t leftPadLen = 2 * maxWMlen;
	uint32_t rightPadLen = maxWMlen;
	// by definition the first position with data will be leftPadLen;
	_firstDatpos = leftPadLen;
	// find maximum position within the fragment and assign lastDatpos
	auto max_it = std::max_element(fragPosVec.begin(), fragPosVec.end());
	_lastDatpos = leftPadLen + (*max_it) - 1; // subtract 1 to make it 0-based

	// define the _size taking into account extensions
	_size = _lastDatpos + rightPadLen + 1;

	// fill with -1.0 as NA sentinel; positions with data are overwritten below
	_modProbVec = std::vector<double>(_size, -1.0);
	for(int i = 0; i < (int)fragPosVec.size(); ++i){
		_modProbVec[_firstDatpos + fragPosVec[i] - 1] = modProbVec[i]; // subtract 1 to make it 0-based
	}
}

// get functions
uint32_t fragProtectData::Size() const{
	return _size;
}
uint32_t fragProtectData::Name() const{
	return _fragID;
}


const double fragProtectData::operator [](uint32_t i) const{
	return _modProbVec[i];
}

fragProtectData & fragProtectData::operator = (const fragProtectData & other){
	if (this != &other){
		_fragID = other._fragID;
		_size = other.Size();

		_firstDatpos = other._firstDatpos;
		_lastDatpos = other._lastDatpos;

		_modProbVec = other._modProbVec;
	}
	return *this;
}

