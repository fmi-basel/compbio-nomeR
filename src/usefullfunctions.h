#ifndef _usefullfunc_h_
#define _usefullfunc_h_


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <string.h>
#include <vector>
using namespace std;

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
