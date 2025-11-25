#include "ftpConfig-class.hpp"


ftpConfig::~ftpConfig(){

}

ftpConfig::ftpConfig(size_t molLength) : length(molLength), ftpOccupancyVec(molLength){
	totalOccupied = 0;
	isFool = 0;
}


// check if a segment can fit into the current configuration without overlaps
bool ftpConfig::canFitFtp(const ftpSegment& cFtpSegm){
	if(cFtpSegm.ftpPosIndex < 0 || cFtpSegm.ftpPosIndex + cFtpSegm.ftpWidth > length){
		return 0;
	}

	return ftpOccupancyVec.range_sum(cFtpSegm.ftpPosIndex + 1,
                                  cFtpSegm.ftpPosIndex + cFtpSegm.ftpWidth) == 0;

}

// add footprint if there is free place and returns TRUE, otherwise returns FALSE
bool ftpConfig::addFtp(const ftpSegment& cFtpSegm){
	if(cFtpSegm.ftpPosIndex < 0 || cFtpSegm.ftpPosIndex + cFtpSegm.ftpWidth > length){
		return 0;
	}

	// check if cFtpSegm can fit in the current configuration
	if(canFitFtp(cFtpSegm)){
		// register segment in the occupancy vector
		// NOTE: the indexing in the occupancy vector is from [1..length]
		ftpOccupancyVec.range_add(cFtpSegm.ftpPosIndex + 1,
                            cFtpSegm.ftpPosIndex + cFtpSegm.ftpWidth, 1);

		// add details of the segment
		// int64_t ftpPosIndex;
		// int64_t ftpStart;
		// int32_t ftpWidth;
		// string ftpGroup;
		// double ftpGroupStartProb;
		// string ftpName;
		// double ftpNameStartProb;

		ftpStartPosVec.push_back(cFtpSegm.ftpStart);    // starting positions of footprints in a configuration
		ftpWidthVec.push_back(cFtpSegm.ftpWidth);       // widths of footprints
		ftpGroupStartProbVec.push_back(cFtpSegm.ftpGroupStartProb);  // start probability for a footprint group, i.e. aggregate across all footprints in the same group
		ftpNameStartProbVec.push_back(cFtpSegm.ftpNameStartProb);   // start probability for the particular footprint with width W
		ftpNameVec.push_back(cFtpSegm.ftpName);         // footprint names
		ftpGroupVec.push_back(cFtpSegm.ftpGroup);        // footprint groups

		// add into total occupancy and set occupancy state

		totalOccupied += cFtpSegm.ftpWidth;
		if(totalOccupied >= length){
			isFool = 1;
		}
		return 1;
	}
	// otherwise return false
	return 0;
}

void ftpConfig::fillConfigVector(vector<int32_t >& cIntSchedFragPos,
                                 vector<int32_t >& cIntSchedFtpWidth,
                                 vector<string >& cIntSchedFtpName,
                                 vector<string >& cIntSchedFtpGroup,
                                 vector<double >& cIntSchedFtpProb){
	cIntSchedFragPos = ftpStartPosVec;
	cIntSchedFtpWidth = ftpWidthVec;
	cIntSchedFtpName = ftpNameVec;
	cIntSchedFtpGroup = ftpGroupVec;
	cIntSchedFtpProb = ftpGroupStartProbVec;
}
