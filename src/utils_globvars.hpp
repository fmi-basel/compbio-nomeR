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

// enumerator for map types in the container
enum fragMapType{GCH, WCG, bisC, otherC,
                 allC // this is used only to get data matrix with all Cs
                 };

// enumerator for type of mismatches
enum seekMismType{CtoT,GtoA};


// global variables for GCH, WCG and bisC contexts, map: context -> strand

// contexts for M.cviPI methylation with endogeous CpG excluded
static unordered_map<string, short int> GCH_CONTEXTS = {
	// on "+" strand
	{"GCA", 1},
	{"GCC", 1},
	{"GCT", 1},
	// on "-" strand
	{"TGC", -1},
	{"GGC", -1},
	{"AGC", -1}
};
// contexts for endogenous CpG with M.cviPI off-targets excluded

static unordered_map<string, short int> WCG_CONTEXTS = {
	// on "+" strand
	{"ACG", 1},
	{"TCG", 1},
	// on "-" strand
	{"CGT", -1},
	{"GCA", -1}
};


// contexts which we use for inference of bisulfite efficiency
// total 16 combinations of NCN
// excluded for bisC
// a <- c("ACG","CCG","GCG","TCG", ## exclude CpGs
//  "GCA","GCC","GCG","GCT", ## exclude GpCs
//  "CCT","CCC","CCA" ## other comb
//)
// allowed combs:
// bisC_allowed <- c("ACA", "ACC", "ACT", "TCA", "TCC", "TCT")
static unordered_map<string, short int> BISC_CONTEXTS = {
	// on "+" strand
	{"ACA", 1},
	{"ACC", 1},
	{"ACT", 1},
	{"TCA", 1},
	{"TCC", 1},
	{"TCT", 1},
	// on "-" strand
	{"TGT", -1},
	{"GGT", -1},
	{"AGT", -1},
	{"TGA", -1},
	{"GGA", -1},
	{"AGA", -1},
};


// SAM spec:
// The case-insensitive base codes ‘=ACMGRSVTWYHKDBN’ are mapped to [0, 15] respectively with
// all other characters mapping to ‘N’ (value 15).

static unordered_map<int, string> nt16_table = {
	{0,"="},
	{1,"A"},
	{2,"C"},
	{3,"M"},
	{4,"G"},
	{5,"R"},
	{6,"S"},
	{7,"V"},
	{8,"T"},
	{9,"W"},
	{10,"Y"},
	{11,"H"},
	{12,"K"},
	{13,"D"},
	{14,"B"},
	{15,"N"}
};

//


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