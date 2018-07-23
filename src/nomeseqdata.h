#ifndef _fastafile_h_
#define _fastafile_h_

#include "sequence.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <Rcpp.h>

using namespace std;

class NOMeSeqData {
  int _numofseq;
  vector<Sequence > _sequences;
//  string _filename;
  int _totallength;
public:
  NOMeSeqData();
  
  ~NOMeSeqData();
  int Size() const;
  int TotalLength() const;
  
  /*
  NOMeSeqData(const vector<string > seq_strings,
              const vector<string > seq_names);
  */
  
  bool create(Rcpp::List _seq_info,
              int maxWMlen);
  
  
  /*
  void NOMeSeqData(const string flnm,int maxwmlen);
  void ReadFasta(const string flnm,int maxwmlen);
   */
  Sequence & operator[](int index);
  void Add(Sequence &);
  void Add(string name, string seq, int wmmaxlen);
  void PrintNames() const;

  void clear();

};

#endif
