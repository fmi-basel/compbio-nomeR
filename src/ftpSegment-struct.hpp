#ifndef _ftpsegment_hpp_
#define _ftpsegment_hpp_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <time.h>
#include <Rcpp.h>
#include <limits>


using namespace std;

struct ftpSegment {
	int32_t ftpPosIndex;
	int32_t ftpStart;
	int32_t ftpWidth;
	string ftpGroup;
	double ftpGroupStartProb;
	string ftpName;
	double ftpNameStartProb;
};


#endif
