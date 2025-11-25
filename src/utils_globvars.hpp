#ifndef _utils_globvars_
#define _utils_globvars_

#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <Rcpp.h>

using namespace std;

// enumerator for method to determine footprint configurations
enum ftpConfigAlgo{POFP,
                 VITERBI
};




// convert string to upper case
string stringToUpper(string strToConvert);



void mychomp(char *s);
void Tokenize(const string& str,
              vector<string>& tokens,
              const string& delimiters );
bool IsNumber(string text);
char revbase(char base);
//char *get_revcomp_seq(const char *seq);
string get_revcomp_seq(const string &seq);

int theta(int pos);
bool compare_vectors ( vector<double> v1,vector<double> v2);
bool compare_vectors_by_column ( const vector<double> v1,const vector<double> v2);
void mySort(vector<vector<double > > & arr,int index);

int letter2index(char letter);
int letter2index_NOMe(char letter);

#endif
