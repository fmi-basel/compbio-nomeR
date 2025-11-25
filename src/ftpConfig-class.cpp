#include "ftpConfig-class.hpp"


ftpConfig::~ftpConfig(){

}

ftpConfig::ftpConfig(size_t molLength) : length(molLength), ftpOccupancyVec(molLength){
	totalOccupied = 0;
	isFool = 0;
}


// check if a segment can fit into the current configuration without overlaps
bool ftpConfig::canFitFtp(const ftpSegment& cFtpSegm,
                          const DNAbind_obj_vector& ftpModels){

	if(cFtpSegm.ftpPosIndex < 0 || cFtpSegm.ftpPosIndex + ftpModels[cFtpSegm.ftpNameIndex]->len > length){
		return 0;
	}
	return ftpOccupancyVec.range_sum(cFtpSegm.ftpPosIndex + 1,
                                  cFtpSegm.ftpPosIndex + ftpModels[cFtpSegm.ftpNameIndex]->len) == 0;
}

// add footprint if there is free place and returns TRUE, otherwise returns FALSE
bool ftpConfig::addFtp(const ftpSegment& cFtpSegm,
                       const DNAbind_obj_vector& ftpModels,
                       const vector<int32_t >& posVecStartProb,
                       vector<int32_t >& cIntSchedFragPos,
                       vector<int32_t >& cIntSchedFtpWidth,
                       vector<string >& cIntSchedFtpName,
                       vector<string >& cIntSchedFtpGroup,
                       vector<double >& cIntSchedFtpProb){
	if(cFtpSegm.ftpPosIndex < 0 || cFtpSegm.ftpPosIndex + ftpModels[cFtpSegm.ftpNameIndex]->len > length){
		return 0;
	}

	// check if cFtpSegm can fit in the current configuration
	if(canFitFtp(cFtpSegm, ftpModels)){
		// register segment in the occupancy vector
		// NOTE: the indexing in the occupancy vector is from [1..length]
		ftpOccupancyVec.range_add(cFtpSegm.ftpPosIndex + 1,
                            cFtpSegm.ftpPosIndex + ftpModels[cFtpSegm.ftpNameIndex]->len, 1);

		// add details of the segment
		cIntSchedFragPos.push_back(posVecStartProb[cFtpSegm.ftpPosIndex]);     // starting positions of footprints in a configuration
		cIntSchedFtpWidth.push_back(ftpModels[cFtpSegm.ftpNameIndex]->len);    // widths of footprints
		cIntSchedFtpName.push_back(ftpModels[cFtpSegm.ftpNameIndex]->name);    // footprint names
		cIntSchedFtpGroup.push_back(ftpModels[cFtpSegm.ftpNameIndex]->group);  // footprint groups
		cIntSchedFtpProb.push_back(cFtpSegm.ftpGroupStartProb);                // start probability for a footprint group, i.e. aggregate across all footprints in the same group

		// ftpStartPosVec.push_back(cFtpSegm.ftpStart);    // starting positions of footprints in a configuration
		// ftpWidthVec.push_back(cFtpSegm.ftpWidth);       // widths of footprints
		// ftpGroupStartProbVec.push_back(cFtpSegm.ftpGroupStartProb);  // start probability for a footprint group, i.e. aggregate across all footprints in the same group
		// ftpNameStartProbVec.push_back(cFtpSegm.ftpNameStartProb);   // start probability for the particular footprint with width W
		// ftpNameVec.push_back(cFtpSegm.ftpName);         // footprint names
		// ftpGroupVec.push_back(cFtpSegm.ftpGroup);        // footprint groups

		// add into total occupancy and set occupancy state

		totalOccupied += ftpModels[cFtpSegm.ftpNameIndex]->len;
		if(totalOccupied >= length){
			isFool = 1;
		}
		return 1;
	}
	// otherwise return false
	return 0;
}

