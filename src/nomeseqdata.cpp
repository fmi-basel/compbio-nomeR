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

bool NOMeSeqData::create(Rcpp::List _seq_list,
                         Rcpp::CharacterVector _seqnames,
                         int maxWMlen){
  _numofseq = 0;
  _totallength = 0;
  
  if(_seq_list.size() != _seqnames.size())
    Rcpp::stop("NOMeSeqData::create: length of sequences list is not equal to length of vector with sequences names.\n");
  
  for(int seq=0; seq < _seq_list.size(); ++seq){
    // get IntegerVector for current sequence
    Rcpp::IntegerVector seq_vec = Rcpp::as<Rcpp::IntegerVector >(_seq_list[seq]);
    // convert IntegerVector to std::vector<unsigned short>
    vector<unsigned short> sq = Rcpp::as<vector<unsigned short> >(seq_vec);
    // get sequence name
    string nm = Rcpp::as<string >(_seqnames[seq]);
    //add sequence into object
    Add(nm,sq,maxWMlen);
    
    }

  return 1;
}

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

void NOMeSeqData::Add(string name, vector<unsigned short> seq, int wmmaxlen)
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
  // however output S will be stating from 1
  // create output vector of vectors
  vector<vector<int > > freq_mat;
  //Rcpp::Rcout<<"Counting frequencies for spacings!"<<endl;
  
  // for each spacing 0:maxSpacing
  for(int s=0; s < maxSpacing; ++s){
    vector<int > freq_vec(10,0); // columns are Spacing, 0,0; 0,1; 0,2; 1,0; 1,1; 1,2; 2,0; 2,1; 2,2
    freq_vec[0] = s + 1;
    // for each sequence
    for(int seq=0; seq < _numofseq; ++seq){
      Sequence sequence = _sequences[seq];
      int seq_len = sequence.Size();
      // go from first position to the last - s + 1
      // note that we go across all position, i.e. including those that we artifically added and filled with NAs(2's)
      for(int pos=0; pos < seq_len - s; pos++){
        //Rcpp::Rcout<<"Spacing= "<<s<<"; seq= "<<seq<<"; pos="<<pos<<endl;
        int letter_pos = sequence[pos];
        int letter_spac = sequence[pos + s];
        //Rcpp::Rcout<<"Spacing= "<<s<<"; seq= "<<seq<<"; pos="<<pos<<"; let_pos="<<letter_pos<<"; let_spac="<<letter_spac<<endl;
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
        } else if(letter_pos == 2 && letter_spac == 0){
          freq_vec[7]++;
        } else if(letter_pos == 2 && letter_spac == 1){
          freq_vec[8]++;
        }
        else if(letter_pos == 2 && letter_spac == 2){
          freq_vec[9]++;
        }
        else {
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

