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
      Rcpp::Rcerr<<"NOMeSeqData::create: Error! At least one element in list of sequences does not contain element DATA\n";
      return(0);
    }
    if(!seqinf.containsElementNamed("NAME")){
      Rcpp::Rcerr<<"NOMeSeqData::create: Error! At least one element in list of sequences does not contain element NAME\n";
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
	Rcpp::Rcout<<"Name: "<<_sequences[i].Name()<<"\tSize: "<<_sequences[i].Size()<<endl;
  }
  Rcpp::Rcout<<"Total length: "<<_totallength<<endl;
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

vector<vector<int> > NOMeSeqData::count_freq_for_spacings(int maxSpacing) const
{
  // here the spacing 0 means that positions are adjacent and gap between them is 0
  // create output vector of vectors
  vector<vector<int> > freq_mat;
  Rcpp::Rcout<<"Counting frequencies for spacings!"<<endl;
  // for each spacing 0:maxSpacing
  for(int s=0; s<=maxSpacing; ++s){
    vector<int > freq_vec(7,0); // columns are Spacing, 0,0; 0,1; 0,2; 1,0; 1,1; 1,2.
    freq_vec[0] = s;
    // for each sequence
    for(int seq=0; seq < _numofseq; ++seq){
      Sequence sequence = _sequences[seq];
      int seq_len = sequence.Size();
      // go from first position to the last - s + 1
      // note that we go across all position, i.e. including those that we artifically added and filled with NAs(2's)
      for(int pos=0; pos < seq_len - s - 10; pos++){
        Rcpp::Rcout<<"Spacing= "<<s<<"; seq= "<<seq<<"; pos="<<pos<<endl;
        int letter_pos = sequence[pos];
        int letter_spac = sequence[pos + s - 1];
        Rcpp::Rcout<<"Spacing= "<<s<<"; seq= "<<seq<<"; pos="<<pos<<"; let_pos="<<letter_pos<<"; let_spac="<<letter_spac<<endl;
        // int letter_pos = 0;
        // int letter_spac = 0;

        if(letter_pos == 0 && letter_spac == 0){
          freq_vec[1]++;
        } else if(letter_pos == 0 && letter_spac == 1){
          freq_vec[2]++;
        } else if(letter_pos == 0 && letter_spac == 2){
          freq_vec[3]++;
        } else if(letter_pos == 1 && letter_spac == 0){
          freq_vec[4]++;
        } else if(letter_pos == 1 && letter_spac == 1){
          freq_vec[5]++;
        } else if(letter_pos == 1 && letter_spac == 2){
          freq_vec[6]++;
        } else {
          
          
          Rcpp::Rcerr<<"NOMeSeqData::count_freq_for_spacings: ERROR! Unknown combination "<<letter_pos<<","<<letter_spac<<" encountered!\n";
          Rcpp::stop("");
          // 
          // exit(1);
        }

      }
    }

    freq_mat.push_back(freq_vec);
  }
  
  return freq_mat;
}

Rcpp::List NOMeSeqData::R_export_spacing_freq(int maxSpacing) const
{
  
  Rcpp::Rcout<<"In R_export_spacing_freq"<<endl;
  vector<vector<int> > freq_mat = count_freq_for_spacings(maxSpacing);
  // 
  // vector<int > spacings;
  // vector<int > freq00;
  // vector<int > freq01;
  // vector<int > freq02;
  // vector<int > freq10;
  // vector<int > freq11;
  // vector<int > freq12;
  // 
  // for(int i=0; i<freq_mat.size(); ++i){
  //   spacings.push_back(freq_mat[i][0]);
  //   freq00.push_back(freq_mat[i][1]);
  //   freq01.push_back(freq_mat[i][2]);
  //   freq02.push_back(freq_mat[i][3]);
  //   freq10.push_back(freq_mat[i][4]);
  //   freq11.push_back(freq_mat[i][5]);
  //   freq12.push_back(freq_mat[i][6]);
  // }
  // 
  // Rcpp::List export_list;
  // 
  // export_list.push_back(Rcpp::wrap(spacings),"Spacing");
  // export_list.push_back(Rcpp::wrap(freq00),"N(0,0)");
  // export_list.push_back(Rcpp::wrap(freq01),"N(0,1)");
  // export_list.push_back(Rcpp::wrap(freq02),"N(0,NA)");
  // export_list.push_back(Rcpp::wrap(freq10),"N(1,0)");
  // export_list.push_back(Rcpp::wrap(freq11),"N(1,1)");
  // export_list.push_back(Rcpp::wrap(freq12),"N(1,NA)");
  // 
  // return export_list;
  
  Rcpp::List export_list;
  return export_list;
  
}

