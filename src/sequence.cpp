#include "sequence.h"


Sequence::Sequence(){

}

Sequence::Sequence(const Sequence & s){
	name = s.Name();
	size = s.Size();
	
	firstDatpos = s.firstDatpos;
	lastDatpos = s.lastDatpos;
	
	
	_seq = s._seq;
}


Sequence::Sequence(string nm,
                   string sequence,
                   int maxWMlen){
	name = nm;
  
  // add flanking regions to line full off 2's (NA's) to take into account that proteins can start outside amplicon 
  string tmp_str(maxWMlen,'2');
  
  string ext_seq = tmp_str + sequence + tmp_str;
  firstDatpos = maxWMlen;
  lastDatpos = maxWMlen + sequence.length() - 1;
  
	size = ext_seq.length();
	if(_seq.size()>0){
		_seq.clear();
	}
	for(int pos=0;pos<size;++pos){
		_seq.push_back(letter2index_NOMe(ext_seq[pos]));
	}
	
	
	
}

Sequence::~Sequence(){

}

string Sequence::Name() const{
	return name;
}

int Sequence::Size() const{
	return size;	
}

unsigned short Sequence::operator [](int i){
	if(i<0 || i>size-1){
	  
	  Rcpp::Rcerr<<"Sequence::operator[]: The index is out of range: index:"<<i<<" size:"<<size<<endl;
	  Rcpp::stop("");
		// cerr<<"Sequence::operator[]: The index is out of range: index:"<<i<<" size:"<<size<<endl;
		// exit(1);
	}
	return _seq[i];

}

vector<unsigned short> Sequence::subseq(int start,int end){
	if(start<0 || end>size-1){
	 
		Rcpp::Rcerr<<"Sequence::subseq: ERROR! The start or the end is out of range: start:"<<start<<"\tend:"<<end<<"\tsize:"<<size<<endl;
		Rcpp::stop("");
		// exit(1);
	}
	vector<unsigned short> subsq;
	for(int i=start;i<=end;++i){
		subsq.push_back(_seq[i]);
	}
	return subsq;
}

Sequence & Sequence::operator = (const Sequence & other){
	if (this != &other){
		name = other.Name();
		size = other.Size();
		if(_seq.size()>0){
			_seq.clear();
		}
		_seq = other._seq;
		
		firstDatpos = other.firstDatpos;
		lastDatpos = other.lastDatpos;
		
		
	}
	return *this;
}

vector<unsigned short> Sequence::getSequence(){
	return _seq;
}

void Sequence::Print() const{
	for(int i=0;i<size;++i){
		Rcpp::Rcout<<_seq[i];
	}
	Rcpp::Rcout<<endl;
}

