#ifndef _sequence_h_
#define _sequence_h_

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <Rcpp.h>
#include "usefullfunctions.h"

using namespace std;


class Sequence {
  public:
	int size;
	string name;
	vector<unsigned short> _seq;
	
	int firstDatpos;
	int lastDatpos;
	
	
	Sequence();
	Sequence(const Sequence & s);
	Sequence(string nm,
          string sequence,
          int maxWMlen);
	~Sequence();
	string Name() const;
	int Size() const;
	unsigned short operator [](int i);
	vector<unsigned short> subseq(int start,int end);
	Sequence & operator = (const Sequence & other);
	vector<unsigned short> getSequence();
	void Print() const;
};

#endif

