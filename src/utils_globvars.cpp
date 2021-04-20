#include "utils_globvars.hpp"




string stringToUpper(string strToConvert)
{
	std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);
	return strToConvert;
}


int _COLUMN_IND_ = 0;


bool compare_vectors_by_column ( const vector<double> v1,const vector<double> v2) {
  extern int  _COLUMN_IND_;
  //Rcpp::Rcout<<"Index for sorting "<<_COLUMN_IND_<<endl;
  if(v1.size() != v2.size()){
    
    Rcpp::stop("compare_vectors_acc_correlation: Vectors must have the same size\n");
    // cerr<<"compare_vectors_acc_correlation: Vectors must have the same size\n";
    // exit(1);
  }
  if( _COLUMN_IND_ == -1){	// means cort by last column
    _COLUMN_IND_ = v1.size() - 1;
  }
  if(  _COLUMN_IND_> v1.size()-1 ||  _COLUMN_IND_ < 0){
    Rcpp::Rcerr<<"Column index is out of range: "<<_COLUMN_IND_<<endl;
    Rcpp::stop("");
    // exit(1);
  }
  if(v1[_COLUMN_IND_] == v2[_COLUMN_IND_]){
    if(_COLUMN_IND_ == 0){
      for(int i=0; i<v1.size();++i){
        if(v1[i] != v2[i]){
          return v1[i] < v2[i];
        }
      }
    }
    else if(_COLUMN_IND_ == v1.size() - 1){
      for(int i=v1.size() - 1; i>0;i--){
        if(v1[i] != v2[i]){
          return v1[i] < v2[i];
        }
      }
    }
  }
  
  return (v1[_COLUMN_IND_] < v2[_COLUMN_IND_]);
}

void mySort(vector<vector<double > > &arr,int index){ // sorts a vector by column index 	ascendently
  extern  int _COLUMN_IND_;
  _COLUMN_IND_ = index;
  //stable_sort(arr.begin(),arr.end(),compare_vectors_by_column);
  sort(arr.begin(),arr.end(),compare_vectors_by_column);
}





void mychomp(char *s) {
  if (s[strlen(s)-1]=='\n')
    s[strlen(s)-1]='\0';
  if (s[strlen(s)-1]=='\r')
    s[strlen(s)-1]='\0';
  
  return;
};

void Tokenize(const string& str,
              vector<string>& tokens,
              const string& delimiters )
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (string::npos != pos || string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
};

bool IsNumber(string text)
{
  if(text.empty()) return false;
  string::const_iterator pos = text.begin();
  while(pos < text.end())
    
    
  {
    if(!isdigit(*pos)) return false;
    ++pos;
  }
  return true;
};

string get_revcomp_seq(const string &seq){
  int pos=0;
  int length=seq.size();
  string revcomp_seq(length,'\0');
  //revcomp_seq = (char *) malloc (length);
  int revpos=0;
  for(pos=length-1;pos>=0;pos--){
    revcomp_seq[revpos] = revbase(seq[pos]);
    revpos++;
  }
  //revcomp_seq[length]='\0';
  return revcomp_seq;
};

/*string get_revcomp_seq(string &seq){
 string revcomp(get_revcomp_seq(seq.c_str()));
 return revcomp;
}
 */
char revbase(char base) {
  if (base == 'A')
    return 'T';
  else if (base == 'C')
    return 'G';
  else if (base == 'G')
    return 'C';
  else if (base == 'T')
    return 'A';
  else if (base == 'N')
    return 'N';
  else if (base == 'W')
    return 'W';
  else if (base == 'S')
    return 'S';
  else if (base == 'Y')
    return 'R';
  else if (base == 'R')
    return 'Y';
  else if (base == 'X')
    return 'X';
  else if (base == 'a')
    return 't';
  else if (base == 'c')
    return 'g';
  else if (base == 'g')
    return 'c';
  else if (base == 't')
    return 'a';
  else if (base == 'n')
    return 'n';
  else if (base == 'w')
    return 'w';
  else if (base == 's')
    return 's';
  else if (base == 'y')
    return 'r';
  else if (base == 'r')
    return 'y';
  else if (base == 'x')
    return 'x';
  else
    return 'N';
};

int theta(int pos){
  int answer;
  if(pos>0){
    answer = 1;
  }
  else{
    answer = 0;
  }
  return answer;
};


//Function for sorting 2D vectors
bool compare_vectors ( vector<double> v1,vector<double> v2) {
  //return (v1[motevo_param.INDEX_OF_WM]<v2[motevo_param.INDEX_OF_WM]);
  return (v1[0]<v2[0]);
  
};

int letter2index(char letter){
  letter = toupper(letter);
  int position;
  if(letter == 'A'){
    position = 0;
  }
  else if(letter == 'C'){
    position = 1;
  }
  else if(letter == 'G'){
    position = 2;
  }
  else if(letter == 'T'){
    position = 3;
  }
  else if(letter == 'X'){
    position = 4;
    Rcpp::Rcout<<"WARNING: The sequence consists letters X\n";
  }
  else{
    Rcpp::Rcerr<<"letter2index:Error! Unknown letter in the sequence "<< letter <<"\n";
    Rcpp::stop("");
    // exit(1);
  }
  return position;
}


int letter2index_NOMe(char letter){
  letter = toupper(letter);
  int position;
  if(letter == '0'){
    position = 0;
  }
  else if(letter == '1'){
    position = 1;
  }
  else if(letter == '2'){
    position = 2;
  }
  else{
    
    
    Rcpp::Rcerr<<"letter2index_NOMe:Error! Unknown letter in the sequence "<< letter <<". Only 0,1,2 allowed!\n";
    Rcpp::stop("");
    // exit(1);
  }
  return position;
}