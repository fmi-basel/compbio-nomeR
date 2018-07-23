#include "nomeseqdata.h"

NOMeSeqData::NOMeSeqData()
{
  _numofseq = 0;
//  _filename = '\0';
  _totallength = 0;

}


bool NOMeSeqData::create(Rcpp::List _seq_info,
                         int maxWMlen){
  _numofseq = 0;
  _totallength = 0;
  
  for(int seq=0; seq < _seq_info.size();++seq){
    Rcpp::List seqinf = Rcpp::as<Rcpp::List >(_seq_info[seq]);
    if(!seqinf.containsElementNamed("DATA")){
      cerr<<"NOMeSeqData::create: Error! At least one element in list of sequences does not contain element DATA\n";
      return(0);
    }
    if(!seqinf.containsElementNamed("NAME")){
      cerr<<"NOMeSeqData::create: Error! At least one element in list of sequences does not contain element NAME\n";
      return(0);
    }
    
    string sq = seqinf["DATA"];
    string nm = seqinf["NAME"];
    
    Add(nm,sq,maxWMlen);
  }
  return 1;
}


/*
FastaFile::FastaFile(const string flnm, int maxwmlen)
{
  ReadFasta(flnm,maxwmlen);
}




void FastaFile::ReadFasta(const string flnm,int maxwmlen)
{
  _filename = flnm;
  ifstream inFile(flnm.c_str());
  if ( !inFile ) {
    cerr << "Can't read file " << flnm <<"\n";
    exit(1);
  }
  string line;
  _totallength = 0;
  while(getline( inFile, line, '\n' )){
	if(line[0] == '#'){ // skip the comments
		continue;
	}
	if(line[0] == '>'){ //the name of the sequence;
		
		string name = line.substr(1,line.size()-1);
		//start searching the nearest sequence
		line.clear();
		while(line.size()==0){ //find first nonempty string
			if(!getline( inFile, line, '\n' )){ // if we reach the end of file
				cerr<<"Couldn't find a sequence for "<<name<<endl;
				exit(1);
			}
		}

    // add flanking regions to line full off 2's (NA's) to take into account that proteins can start outside amplicon 
    string tmp_str(maxwmlen,'2');

		Sequence newseq(name,tmp_str + line + tmp_str);
		_totallength += newseq.Size();
		_sequences.push_back(newseq);
	}
  }
  _numofseq = _sequences.size();
  firstDatPos = maxwmlen;
}
*/



NOMeSeqData::~NOMeSeqData()
{
  
}

int NOMeSeqData::Size() const
{
  return _numofseq;
}

int NOMeSeqData::TotalLength() const
{
  return _totallength;
}


void NOMeSeqData::PrintNames() const
{
  for(int i=0;i<_numofseq;++i){
	cout<<"Name: "<<_sequences[i].Name()<<"\tSize: "<<_sequences[i].Size()<<endl;
  }
  cout<<"Total length: "<<_totallength<<endl;
}


Sequence & NOMeSeqData::operator[] (int index)
{
  return _sequences[index];
}

void NOMeSeqData::Add(Sequence &seq)
{
  Sequence newseq(seq);
  _totallength += newseq.Size();
  _sequences.push_back(newseq);
}

void NOMeSeqData::Add(string name, string seq, int wmmaxlen)
{
  Sequence newseq(name,seq,wmmaxlen);
  _totallength += newseq.Size();
  _sequences.push_back(newseq);
  _numofseq++;
}

void NOMeSeqData::clear(){
  _numofseq =0 ;
  _sequences.clear();
  //_sequences.swap(vector<Sequence >());
  //  string _filename;
  _totallength=0;
  
}