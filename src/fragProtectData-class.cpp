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
	_protectVec = s._protectVec;
}

fragProtectData::fragProtectData(const uint32_t fragID,
                                 const vector<uint32_t>& fragPosVec,
                                 const vector<uint8_t>& protectVec,
                                 int maxWMlen){
	_fragID = fragID;
	// here, we have to redefine positions within the fragment to take into account padding by NAs of size maxWMlen
	// by definition the first position with data will be maxWMlen;
	_firstDatpos = maxWMlen;
	// find maximum position within the fragment and assign lastDatpos
	auto max_it = std::max_element(fragPosVec.begin(), fragPosVec.end());
	_lastDatpos = maxWMlen + (*max_it);

	// define the _size taking into account extensions
	_size = _lastDatpos + maxWMlen + 1;

	// add protection data to _protectVec
	_protectVec = std::vector<uint8_t>(_size, 2);
	for(int i = 0; i < fragPosVec.size(); ++i){
		_protectVec[_firstDatpos + fragPosVec[i]] = protectVec[i];
	}
}

// get functions
uint32_t fragProtectData::Size() const{
	return _size;
}
uint32_t fragProtectData::Name() const{
	return _fragID;
}


const uint8_t fragProtectData::operator [](uint32_t i) const{
	return _protectVec[i];
}

fragProtectData & fragProtectData::operator = (const fragProtectData & other){
	if (this != &other){
		_fragID = other._fragID;
		_size = other.Size();

		_firstDatpos = other._firstDatpos;
		_lastDatpos = other._lastDatpos;
		
		_protectVec = other._protectVec;
	}
	return *this;
}

