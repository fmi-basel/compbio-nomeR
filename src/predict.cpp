#include "predict.hpp"

Predict::Predict(){}
Predict::~Predict(){
  // object contains vector of pointers. need to call explicitly clear to delete
  BINDING_OBJECTS.clear();
  }

Predict::Predict(const Rcpp::List& data,
                 const Rcpp::CharacterVector& fragnames,
                 const Rcpp::List& binding_models,
                 const Rcpp::NumericVector& bgprotectprob,
                 const Rcpp::NumericVector& bgprior)
{
  Create(data,
         fragnames,
         binding_models,
         bgprotectprob,
         bgprior);
}

bool Predict::Create(const Rcpp::List& data,
                     const Rcpp::CharacterVector& fragnames,
                     const Rcpp::List& binding_models,
                     const Rcpp::NumericVector& bgprotectprob,
                     const Rcpp::NumericVector& bgprior)
{
  // verbosity
  extern bool _VERBOSE_;
  
  // create object with parameters
  double bgcoverprob_ = Rcpp::as<double >(bgprotectprob);
  double bgprior_ = Rcpp::as<double >(bgprior);
  if(_VERBOSE_){
    Rcpp::Rcout<<"Creating PARAMS object..."<<endl;
  }
  PARAMS.setParams(bgcoverprob_,
                   bgprior_);
  
  
  // create object with background/footprint models
  if(_VERBOSE_){
    Rcpp::Rcout<<"Creating BINDING_OBJECTS object..."<<endl;
  }
  BINDING_OBJECTS.create(binding_models,PARAMS);
  
  // create object with NOMeSeq data
  if(_VERBOSE_){
    Rcpp::Rcout<<"Creating SEQUENCES object..."<<endl;
  }
  SEQUENCES.create(data,
                   fragnames,
                   BINDING_OBJECTS.maxwmlen);
  
  
  // Initialize print_indexes
  if(PARAMS.printoutonly != "All"){
    Tokenize(PARAMS.printoutonly,print_names," \t,;{}");
    
    for(int i=0;i<print_names.size();++i){
      
      for(int wm=0;wm<BINDING_OBJECTS.Size();++wm){
        if(BINDING_OBJECTS[wm] -> name == print_names[i]){
          print_indexes.push_back(wm);
        }
      }
    }
    if(print_indexes.size() == 0){
      Rcpp::stop("Forward_Backward_algorithm:: ERROR! Incorrect value of parameter printoutonly. No objects to calculate were found.\n");
    }
  }
  else{
    for(int wm=0;wm< BINDING_OBJECTS.Size();++wm){
      print_indexes.push_back(wm);
      vector<string >::iterator it;
      it = find(print_names.begin(),print_names.end(),BINDING_OBJECTS[wm]->name);
      if(it == print_names.end())
        print_names.push_back(BINDING_OBJECTS[wm]->name);
    }
  }
  
  // define maps indexes2names and names2indexes
  // names2indexes - this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  for(int i=0;i<print_names.size();++i){
    vector<int > tmp;
    for(int wm=0;wm<print_indexes.size();++wm){
      if(BINDING_OBJECTS[print_indexes[wm]]->name == print_names[i])
        tmp.push_back(print_indexes[wm]);
    }
    names2indexes.push_back(tmp);
  }
  
  // define names2indicesinprobarray
  for(int i=0;i<print_names.size();++i){
    vector<int > tmp;
    for(int wm=0;wm<print_indexes.size();++wm){
      if(BINDING_OBJECTS[print_indexes[wm]]->name == print_names[i])
        tmp.push_back(wm);
    }
    names2indicesinprobarray.push_back(tmp);
    
  }
  
  // Allocate memory for F, R and Prob
  int maxwmlen = BINDING_OBJECTS.maxwmlen;
  for(int seq=0; seq < SEQUENCES.Size();++seq){
    vector<double > tmp(SEQUENCES[seq].Size() + maxwmlen + 1,1);
    F.push_back(tmp);
    tmp.resize(SEQUENCES[seq].Size() + 2,1);
    R.push_back(tmp);
    vector<double > tmp2(SEQUENCES[seq].Size() + 1, 0);
    vector<vector<double > > tmp1(print_indexes.size(),tmp2);
    Prob.push_back(tmp1);
    
  }
  // Allocate genomesummary
  vector<double > tmp(7,0);
  genomesummary.resize(print_names.size(),tmp);
  return 1;
}


/* // [[Rcpp::depends(RcppProgress)]]
 * */
bool Predict::Run(int ncpu) 
{
  extern bool _VERBOSE_;
  
  int numberofobjects = BINDING_OBJECTS.Size();
  int seq = 0;
  
  // initial value for partition sums.
  // when footprint priors are normalized, i.e. sum of all priors is 1, then initial values for partition sums is always 1.
  // we normalize the priors, therefore we set value to 1.
  double part_init = 1;
  
#ifdef _OPENMP
  omp_set_nested(true);
  omp_set_num_threads(ncpu);
  if(_VERBOSE_)
    Rcpp::Rcout<<"Running prediction with "<<omp_get_max_threads()<<" cpu."<<endl;
#endif
  
#pragma omp parallel private(seq)
{
  
#pragma omp for schedule(dynamic)
  for(seq=0;seq < SEQUENCES.Size();++seq){
    
    // check for user interruption every 10 sequences
    //if(seq % 10 == 0){
    // if(Progress::check_abort()){
    //   //Rcpp::Rcout<<"Interupted by user"<<endl;
    //   return(0);
    // }
    //}
    
    
    int seqlength = SEQUENCES[seq].Size();
    // calculate forward partition summ
    F[seq][0] = part_init;
    vector<double > pf(numberofobjects,1);
    for(int pos=1;pos <= seqlength + BINDING_OBJECTS.maxwmlen; ++pos){
      double summ=0;
      for(int wm=0;wm < numberofobjects; ++wm){
        pf[wm] = 1;
        int objlen = BINDING_OBJECTS[wm]->len;
        
        if(BINDING_OBJECTS[wm]->prior > 0){
          pf[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq, pos - objlen);
          for(int i = pos - objlen + 1; i <= pos - 1; ++i){
            
            if(i>=0){
              pf[wm] *= F[seq][i];
            } else {
              pf[wm] *= part_init;
            }
          }
          
        }
        else{
          pf[wm] = 0;
        }
        summ += pf[wm];
      }
      
      F[seq][pos] = 1/summ;
      if(pos<=seqlength){
        for(int wm=0;wm<print_indexes.size();++wm){
          Prob[seq][wm][pos] = pf[print_indexes[wm]];
        }
      }
    }
    
    // calculate backward partition summ
    vector<double > pb(numberofobjects, 1);
    R[seq][seqlength + 1] = part_init;
    for(int pos = seqlength; pos>=1;--pos){
      double summ = 0;
      for(int wm = 0;wm < numberofobjects; ++wm){
        pb[wm] = 1;
        int objlen = BINDING_OBJECTS[wm]->len;
        if(BINDING_OBJECTS[wm]->prior > 0){
          pb[wm] = BINDING_OBJECTS[wm]->get_score(SEQUENCES,seq,pos - 1);
          for(int i = pos+1; i <= pos+objlen-1; ++i){
            if(i<=seqlength + 1){
              pb[wm] *= R[seq][i];
            } else{
              pb[wm] *= part_init;
            }
          }
          
        }
        else{
          pb[wm] = 0;
        }
        summ += pb[wm];
      }
      
      R[seq][pos] = 1/summ;
    }
    // сalculate Z = Fn/Rn
    
    // calculate initial value for Z based on requirement that total coverage at L must be 1
    double zsumm =0;
    for(int wm = 0;wm < numberofobjects; ++wm){
      int objlen = BINDING_OBJECTS[wm]->len;
      double wmsumm = 0;
      for(int pos=seqlength; pos<= seqlength + objlen - 1; ++pos){
        double prod = 1;
        for(int j=pos - objlen + 1; j<=seqlength;++j){
          prod *= F[seq][j];
        }
        wmsumm += pow(part_init,pos - seqlength) * prod;
      }
      zsumm += BINDING_OBJECTS[wm]->prior * wmsumm;
    }
    
    double z_init = 1/zsumm;
    
    R[seq][seqlength + 1] = z_init;
    for(int pos = seqlength;pos >= 1;--pos){
      R[seq][pos] = F[seq][pos] * R[seq][pos + 1]/R[seq][pos];
    }
    // calculate posteriors
    for(int pos = 1;pos <= seqlength;++pos){
      
      for(int wm=0;wm < print_indexes.size(); ++wm){
        int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
        if(pos + objlen - 1 <= seqlength)
          Prob[seq][wm][pos] = Prob[seq][wm][pos + objlen - 1] * F[seq][pos + objlen -1] * R[seq][pos+objlen];
        else
          Prob[seq][wm][pos] = 0;
      }
    }
  }
  
}

return(1);
}


vector<double > Predict::getPriors()
{
  int numberofobjects = BINDING_OBJECTS.Size();
  vector<double > curpriors(numberofobjects,0);
  
  double totsum=0;
  for(int wm=0;wm<numberofobjects;++wm){
    for(int seq =0;seq<Prob.size();++seq){
      for(int pos = 1;pos <= SEQUENCES[seq].Size();++pos){
        curpriors[wm] +=Prob[seq][wm][pos];
        totsum += Prob[seq][wm][pos];
      }
    }
  }
  
  if(totsum == 0){
    Rcpp::stop("Predict::getPriors(): Error! total sum of probabilities is 0! Exit.\n");
  }
  for(int wm=0;wm<numberofobjects;++wm){
    curpriors[wm] /= totsum;
  }
  
  return(curpriors);
}


void Predict::SetGenomeSummary(){
  
  for(int i=0;i<genomesummary.size();++i){
    for(int j=0;j<genomesummary[i].size();++j)
      genomesummary[i][j] = 0;
  }

  for(int name=0;name<print_names.size();++name){
    {
      
      for(int seq =0;seq<Prob.size();++seq){
        
        for(int pos = 1;pos <= SEQUENCES[seq].Size();++pos){
          double totalprob=0;
          for(int wm =0;wm<names2indicesinprobarray[name].size();++wm)
            totalprob += Prob[seq][names2indicesinprobarray[name][wm]][pos];
          
          genomesummary[name][0] += totalprob;
          genomesummary[name][1] += totalprob;
          genomesummary[name][2] += BINDING_OBJECTS[names2indexes[name][0]]->len * totalprob/SEQUENCES.TotalLength();
          if(totalprob < 0.5){
            genomesummary[name][3] += totalprob;
            genomesummary[name][5]++;
          }
          else{
            genomesummary[name][4] += totalprob;
            genomesummary[name][6]++;
          }
        }
      }
      
    }
  }
  // calculate priors as ratio between number of sites and total number of sites for all binding objects
  double allobjtot = 0;
  for(int wm=0; wm<BINDING_OBJECTS.Size();++wm){
    for(int seq =0;seq<Prob.size();++seq){
      for(int pos = 1;pos <= SEQUENCES[seq].Size();++pos){
        allobjtot += Prob[seq][wm][pos];
      }
    }
  }
  
  for(int name=0;name<print_names.size();++name){
    genomesummary[name][0] /= allobjtot;
  }
  
}



double Predict::get_coverage_prob_at_pos_name_index(int seq, int pos, int name_index){ //wm is the index in print_names
  
  if(Prob.empty()){
    Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
  }

  if(seq < 0 || seq>SEQUENCES.Size() - 1){
    Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The sequence index is out of range.\n");
  }
  
  if(pos<1 || pos>SEQUENCES[seq].Size()){
    Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The position is out of range.\n");
  }
  
  if(name_index<0 || name_index > print_names.size() -1){
    Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! The name_index is out of range.\n");
  }
  
  string name = print_names[name_index];
  double coverage=0;
  int objlen = BINDING_OBJECTS[names2indexes[name_index][0]]->len;
  for(int i=0;i<names2indicesinprobarray[name_index].size();++i){
    for(int p = max(pos - objlen + 1,1);p <= pos; ++p){
      coverage += Prob[seq][names2indicesinprobarray[name_index][i]][p];
    }
  }
  
  if(coverage<0 || coverage>1+1e-5){
    Rcpp::stop("Predict::get_coverage_prob_at_pos_name_index: ERROR! Incorrect value of coverage probability.\n");
  }
  
  return coverage;
}



double Predict::get_coverage_prob_at_pos(int seq, int pos, int wm){ //wm is the index in Prob
  
  if(Prob.empty()){
    Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
  }

  if(seq < 0 || seq>SEQUENCES.Size() - 1){
    Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The sequence index is out of range.\n");
  }
  
  if(pos<1 || pos>SEQUENCES[seq].Size()){
    Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The position is out of range.\n");
  }
  if(wm<0 || wm > Prob[seq].size() -1){
    Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! The wm index is out of range.\n");
  }
  
  double coverage=0;
  int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
  for(int p = max(pos - objlen + 1,1);p <= pos; ++p)
    coverage += Prob[seq][wm][p];
  
  if(coverage < 0 || coverage > 1+1e-5){
    Rcpp::stop("Predict::get_coverage_prob_at_pos: ERROR! Incorrect value of coverage probability.\n");
  }
  return coverage;
  
}



Rcpp::List Predict::getStartProbDF(){
  extern bool _VERBOSE_;
  
  Rcpp::List out_list;
  vector<string > seqnames;
  vector<int > positions;
  // fill seq and pos
  for(int seq=0;seq < SEQUENCES.Size();seq++){
    int firstDatPos = SEQUENCES[seq].firstDatpos;
    int lastDatPos = SEQUENCES[seq].lastDatpos;

    for(int position = firstDatPos;position <= lastDatPos;++position){
      seqnames.push_back(SEQUENCES[seq].Name());
      positions.push_back(position - firstDatPos + 1);
    }
  }
  
  out_list.push_back(Rcpp::wrap(seqnames),"seq");
  out_list.push_back(Rcpp::wrap(positions),"pos");

  // fill start probabilities
  for(int i=0;i<print_names.size();++i){
    vector<double > stprob;
    int seq=0;

    for(seq=0;seq < SEQUENCES.Size();seq++){
      int firstDatPos = SEQUENCES[seq].firstDatpos;
      int lastDatPos = SEQUENCES[seq].lastDatpos;
      
      for(int position = firstDatPos;position <= lastDatPos;++position){
        double totalprob=0;
        for(int j=0;j<names2indicesinprobarray[i].size();++j){
          totalprob += Prob[seq][names2indicesinprobarray[i][j]][position];
        }
        stprob.push_back(totalprob);
      }
    }
    
    out_list.push_back(Rcpp::wrap(stprob),print_names[i]);
  }
  
  return out_list;
}



Rcpp::List Predict::getCoverProbDF(){

  Rcpp::List out_list;
  vector<string > seqnames;
  vector<int > positions;
  // fill seq and pos
  for(int seq=0;seq < SEQUENCES.Size();seq++){
    int firstDatPos = SEQUENCES[seq].firstDatpos;
    int lastDatPos = SEQUENCES[seq].lastDatpos;
    
    for(int position = firstDatPos;position <= lastDatPos;++position){
      seqnames.push_back(SEQUENCES[seq].Name());
      positions.push_back(position - firstDatPos + 1);
    }
  }
  
  out_list.push_back(Rcpp::wrap(seqnames),"seq");
  out_list.push_back(Rcpp::wrap(positions),"pos");

  // fill coverage probabilities
  for(int i=0;i<print_names.size();++i){
    
    vector<double > covprob;
    for(int seq=0;seq < SEQUENCES.Size();seq++){
      int firstDatPos = SEQUENCES[seq].firstDatpos;
      int lastDatPos = SEQUENCES[seq].lastDatpos;
      
      for(int position = firstDatPos;position <= lastDatPos;++position){
        double coverprob=get_coverage_prob_at_pos_name_index(seq,position,i);
        covprob.push_back(coverprob);
      }
    }
    out_list.push_back(Rcpp::wrap(covprob),print_names[i]);
  }
  
  return out_list;
}



Rcpp::List Predict::getGenomeSummaryDF(){
  
  // calculate statistics across the whole amplicon (all fragments)
  SetGenomeSummary();
  Rcpp::List out_list;
  const char *fnms[] = {"Prior",
                        "Expected number of sites",
                        "Coverage",
                        "Sites<0.5",
                        "Sites>=0.5",
                        "Positions<0.5",
                        "Positions>=0.5"};
  vector<string > field_names(fnms,fnms + 7) ;
  out_list.push_back(Rcpp::wrap(field_names),"Statistics");
  for(int wm=0;wm<print_names.size();++wm){
    out_list.push_back(Rcpp::wrap(genomesummary[wm]),print_names[wm]);
  }

  return out_list;
}


void Predict::clear(){
  
  SEQUENCES.clear();
  PARAMS.clear();
  BINDING_OBJECTS.clear();
  
  F.clear();//swap(vector<vector<double> >());
  R.clear();//(vector<vector<double> >());
  Prob.clear();//(vector<vector<vector<double> > >());
  genomesummary.clear();//(vector<vector<double > > ()); //this vector contains expected prior, expected number of sites and expected coverage for the whole genome and for each TF
  
  print_indexes.clear();//swap(vector<int >());	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
  names2indexes.clear();//swap(vector<vector<int > >()); // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  print_names.clear();//swap(vector<string >()); // this array contain names of the objects that will be printed
  
  names2indicesinprobarray.clear();//swap(vector<vector<int > >()); // this array contains map i - index in names to subarray of indices in Prob array
  
}

