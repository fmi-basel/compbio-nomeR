#include "ftpConfig-class.hpp"


ftpConfig::~ftpConfig(){

}

ftpConfig::ftpConfig(size_t molLength)
{
	_occupancy = std::vector<uint8_t>(molLength, 0);
}


// add footprint if there is free place and returns TRUE, otherwise returns FALSE
bool ftpConfig::addFtp(const int64_t& _ftpStart,
                       const int32_t& _ftpWidth,
                       const double& _ftpGroupStartProb,
                       const double& _ftpNameStartProb,
                       const string& _ftpName,
                       const string& _ftpGroup){

	// find closest occupied position after _ftpStart
	int64_t nextOccupPos = findNextOccupPos(_ftpStart);

	if(nextOccupPos >= _ftpStart + _ftpWidth){ // 0-based indexing
		ftpStartPosVec.push_back(_ftpStart);
		ftpWidthVec.push_back(_ftpWidth);
		ftpGroupStartProbVec.push_back(_ftpGroupStartProb);
		ftpNameStartProbVec.push_back(_ftpNameStartProb);
		ftpNameVec.push_back(_ftpName);
		ftpGroupVec.push_back(_ftpGroup);

		// mark occupancy for the new footprint
		markOccupPositions(_ftpStart,_ftpWidth);
    // return true
		return 1;
	}

  // otherwise return false
	return 0;

}
