#include "forward_backward.h"

Forward_Backward_algorithm::Forward_Backward_algorithm()
{

}



bool Forward_Backward_algorithm::Create()
{
	extern parameters PARAMS;
	extern DNAbind_obj_vector BINDING_OBJECTS;
	extern NOMeSeqData SEQUENCES;
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
			// cerr<<"Forward_Backward_algorithm:: ERROR! Incorrect value of parameter printoutonly. No objects to calculate were found.\n";
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

	//Now define maps indexes2names and names2indexes
	//names2indexes - this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
	for(int i=0;i<print_names.size();++i){
		vector<int > tmp;
		for(int wm=0;wm<print_indexes.size();++wm){
			if(BINDING_OBJECTS[print_indexes[wm]]->name == print_names[i])
				tmp.push_back(print_indexes[wm]);
		}
		names2indexes.push_back(tmp);
	}
	
	// now define names2indicesinprobarray

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


Forward_Backward_algorithm::~Forward_Backward_algorithm(){


}

void Forward_Backward_algorithm::Run(vector<double > priors,
                                     int ncpu) 
{

	extern DNAbind_obj_vector BINDING_OBJECTS;
	extern NOMeSeqData SEQUENCES;
	extern bool _VERBOSE_;
	if(!priors.empty()){
		if(priors.size() != BINDING_OBJECTS.numberofobjects){
		  Rcpp::stop("Forward_Backward_algorithm::Run_Low_Mem: ERROR: The size of priors vector doesn't equal the size of BINDING_OBJECTS vector.\n");
		}
		for(int i=0;i<BINDING_OBJECTS.numberofobjects;++i){
			for(int j=0;j<BINDING_OBJECTS.names2index[i].size();++j){
				double prior = priors[i]/BINDING_OBJECTS.names2index[i].size();
				BINDING_OBJECTS[BINDING_OBJECTS.names2index[i][j]] -> prior = prior;
				//BINDING_OBJECTS[BINDING_OBJECTS.names2index[i][j]] -> betta = bettas[i];
			}
		}
	}
	
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
	  
	  //check for interrupt every 10 sequences
	  // if(seq % 1000 == 0)
	  //   Rcpp::checkUserInterrupt();
	  
		int seqlength = SEQUENCES[seq].Size();
	
//#pragma omp parallel shared(seq)
//{
//#pragma omp sections
//{
	
	//#pragma omp section
	//{

		// Forward partition summ
		F[seq][0] = part_init;
		vector<double > pf(numberofobjects,1);
		for(int pos=1;pos <= seqlength + BINDING_OBJECTS.maxwmlen; ++pos){
			double summ=0;
			for(int wm=0;wm < numberofobjects; ++wm){
				pf[wm] = 1;
				int objlen = BINDING_OBJECTS[wm]->len;
				
				if(BINDING_OBJECTS[wm]->prior > 0){
					//Rcpp::Rcout<<pos<<"\t"<<BINDING_OBJECTS[wm]->name<<endl;
					//int rl = pos - objlen;

					//Rcpp::Rcout<<"Seq="<<SEQUENCES[seq].Name()<<"; pos="<<pos<<"; name="<<BINDING_OBJECTS[wm]->name<<"; score="<< BINDING_OBJECTS[wm]->get_score(seq, pos - objlen)<<endl;
					
					pf[wm] = BINDING_OBJECTS[wm]->get_score(seq, pos - objlen);
					for(int i = pos - objlen + 1; i <= pos - 1; ++i){
						//Rcpp::Rcout<<"i="<<i<<"; F="<<F[seq][i]<<endl;
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
			//double bgposscore = BINDING_OBJECTS[numberofobjects - 1]->get_score(seq, pos);

			F[seq][pos] = 1/summ;
			//Rcpp::Rcout<<"Seq="<<SEQUENCES[seq].Name()<<"; pos=\t"<<pos<<"; sum="<<summ<<"; F=\t"<<F[seq][pos]<<endl;
			if(pos<=seqlength){
				for(int wm=0;wm<print_indexes.size();++wm){
					//Rcpp::Rcout<<"pf "<<pf[print_indexes[i]]<<endl;
					Prob[seq][wm][pos] = pf[print_indexes[wm]];
				}
			}
			//Rcpp::Rcout<<"seq="<<SEQUENCES[seq].Name()<<"; pos="<<pos<<"; Pr="<<Pr[1][pos]<<"\n";
		}
	
	
	//#pragma omp section
	//{
		
		// Backward partition summ
		vector<double > pb(numberofobjects, 1);
		R[seq][seqlength + 1] = part_init;
		for(int pos = seqlength; pos>=1;--pos){
			double summ = 0;
			for(int wm = 0;wm < numberofobjects; ++wm){
				pb[wm] = 1;
				int objlen = BINDING_OBJECTS[wm]->len;
				if(BINDING_OBJECTS[wm]->prior > 0){
					pb[wm] = BINDING_OBJECTS[wm]->get_score(seq,pos - 1);
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
			
			//double bgposscore = BINDING_OBJECTS[numberofobjects - 1]->get_score(seq, pos-1);
			
			R[seq][pos] = 1/summ;
			//Rcpp::Rcout<<"Seq="<<SEQUENCES[seq].Name()<<"; pos=\t"<<pos<<"; sum="<<summ<<"; R=\t"<<R[seq][pos]<<endl;
			
		}
	

	//}
//}
//}
		// Calculate Z = Fn/Rn
		
		/* First we have to calculate initial value for Z from the condition that total coverage at L must be 1
		*/
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
			//Rcpp::Rcout<<"pos\t"<<pos<<"\tZ=\t"<<R[seq][pos]<<endl;
		}
		//Calculate posteriors
		for(int pos = 1;pos <= seqlength;++pos){
			//Rcpp::Rcout<<pos;
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
/*	time(&end);
        Rcpp::Rcout<<"Forward_Backward_algorithm::Run_Low_Mem():It took "<<difftime(end,start)<<endl;
*/
}



void Forward_Backward_algorithm::Run(int ncpu) 
{
	vector<double > priors;
	Run(priors,
      ncpu);
}





vector<double > Forward_Backward_algorithm::getPriors(){

	extern DNAbind_obj_vector BINDING_OBJECTS;
  	extern NOMeSeqData SEQUENCES;

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
	  Rcpp::stop("Forward_Backward_algorithm::GetPriors(): Error! total sum of probabilities is 0! Exit.\n");
	  
		// cerr<<"Forward_Backward_algorithm::GetPriors(): Error! total sum of probabilities is 0! Exit.\n";
		// exit(1);
	}
	for(int wm=0;wm<numberofobjects;++wm){
		curpriors[wm] /= totsum;
	}

	return(curpriors);
}


void Forward_Backward_algorithm::SetGenomeSummary(){
  for(int i=0;i<genomesummary.size();++i){
    for(int j=0;j<genomesummary[i].size();++j)
      genomesummary[i][j] = 0;
  }
  extern DNAbind_obj_vector BINDING_OBJECTS;
  extern NOMeSeqData SEQUENCES;
  
  
  
  
  for(int name=0;name<print_names.size();++name){
    //double a0=0,a1=0,a2=0,a3=0,a4=0,a5=0,a6=0;
    
//#pragma omp parallel
{
//#pragma omp for schedule(dynamic) reduction (+:a0,a1,a2,a3,a4,a5,a6)
  for(int seq =0;seq<Prob.size();++seq){
    //int proc_id = 
    //Rcpp::Rcout<<"Process id: "<<omp_get_thread_num()<<endl;
    for(int pos = 1;pos <= SEQUENCES[seq].Size();++pos){
      double totalprob=0;
      for(int wm =0;wm<names2indicesinprobarray[name].size();++wm)
        totalprob += Prob[seq][names2indicesinprobarray[name][wm]][pos];
      
      genomesummary[name][0] += totalprob;
      genomesummary[name][1] += totalprob;
      genomesummary[name][2] += BINDING_OBJECTS[names2indexes[name][0]]->len * totalprob/SEQUENCES.TotalLength();
      // a0 += totalprob;
      // a1 += totalprob;
      // a2 += BINDING_OBJECTS[names2indexes[name][0]]->len * totalprob/SEQUENCES.TotalLength();
      if(totalprob < 0.5){
        genomesummary[name][3] += totalprob;
        genomesummary[name][5]++;
        // a3 += totalprob;
        // a5++;
      }
      else{
        genomesummary[name][4] += totalprob;
        genomesummary[name][6]++;
        // a4 += totalprob;
        // a6++;
      }
      
      
    }
  }
  
}
// genomesummary[name][0] = a0;
// genomesummary[name][1] = a1;
// genomesummary[name][2] = a2;
// genomesummary[name][3] = a3;
// genomesummary[name][5] = a5;
// genomesummary[name][4] = a4;
// genomesummary[name][6] = a6;

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


double Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index(int seq, int pos, int name_index){ //wm is the index in print_names
	if(Prob.empty()){
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
	  
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n";
		// exit(1);
	}
	extern DNAbind_obj_vector BINDING_OBJECTS;
  	extern NOMeSeqData SEQUENCES;
	if(seq < 0 || seq>SEQUENCES.Size() - 1){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The sequence index is out of range.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The sequence index is out of range "<<seq<<endl;
		// exit(1);
	}
	if(pos<1 || pos>SEQUENCES[seq].Size()){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The position is out of range.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The position is out of range "<<pos<<endl;
		// exit(1);
	}
	if(name_index<0 || name_index > print_names.size() -1){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The name_index is out of range.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! The name_index is out of range "<<name_index<<endl;
		// exit(1);
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
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! Incorrect value of coverage probability.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos_name_index: ERROR! Incorrect value of coverage probability: "<<seq<<"\t"<<pos<<"\t"<<name<<"\t"<<coverage<<endl;
		//exit(1);
	}
	return coverage;

}

double Forward_Backward_algorithm::get_coverage_prob_at_pos(int seq, int pos, int wm){ //wm is the index in Prob
	if(Prob.empty()){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! Unable to calculate coverage probability. The vector Prob is empty\n";
		// exit(1);
	}
	extern DNAbind_obj_vector BINDING_OBJECTS;
  	extern NOMeSeqData SEQUENCES;
	if(seq < 0 || seq>SEQUENCES.Size() - 1){
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The sequence index is out of range.\n");
	  
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The sequence indes is out of range "<<seq<<endl;
		// exit(1);
	}
	if(pos<1 || pos>SEQUENCES[seq].Size()){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The position is out of range.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The position is out of range "<<pos<<endl;
		// exit(1);
	}
	if(wm<0 || wm > Prob[seq].size() -1){
	  
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The wm index is out of range.\n");
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! The wm index is out of range "<<wm<<endl;
		// exit(1);
	}
	
	
	double coverage=0;
	int objlen = BINDING_OBJECTS[print_indexes[wm]]->len;
	for(int p = max(pos - objlen + 1,1);p <= pos; ++p)
		coverage += Prob[seq][wm][p];
		
	if(coverage < 0 || coverage > 1+1e-5){
	  Rcpp::stop("Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! Incorrect value of coverage probability.\n");
	  
		// cerr<<"Forward_Backward_algorithm::get_coverage_prob_at_pos: ERROR! Incorrect value of coverage probability: "<<seq<<"\t"<<pos<<"\t"<<BINDING_OBJECTS[print_indexes[wm]]->name<<"\t"<<coverage<<endl;
		// exit(1);
	}
	return coverage;

}



Rcpp::List  Forward_Backward_algorithm::getStartProbDF(){
  extern parameters PARAMS;
  extern NOMeSeqData SEQUENCES;
  extern DNAbind_obj_vector BINDING_OBJECTS;
  extern bool _VERBOSE_;

  Rcpp::List out_list;
  vector<string > seqnames;
  vector<int > positions;
  // fill seq and pos
  for(int seq=0;seq < SEQUENCES.Size();seq++){
    int firstDatPos = SEQUENCES[seq].firstDatpos;
    int lastDatPos = SEQUENCES[seq].lastDatpos;
//    Rcpp::Rcout<<"FirstPos="<<firstDatPos<<"\tLastPos="<<lastDatPos<<endl;
    for(int position = firstDatPos;position <= lastDatPos;++position){
      seqnames.push_back(SEQUENCES[seq].Name());
      positions.push_back(position - firstDatPos + 1);
    }
    
    // for(int position = 1;position <= SEQUENCES[seq].Size();++position){
    //   seqnames.push_back(SEQUENCES[seq].Name());
    //   positions.push_back(position - firstDatPos + 1);
    // }
    
  }
  
  //out_list[0] = Rcpp::wrap(seqnames);
  //out_list[1] = Rcpp::wrap(positions);
  out_list.push_back(Rcpp::wrap(seqnames),"seq");
  out_list.push_back(Rcpp::wrap(positions),"pos");
  //out_list.push_back(Rcpp::as<Rcpp::IntegerVector>(positions),"pos");
  
  // fill start probabilities
  for(int i=0;i<print_names.size();++i){
    //Rcpp::NumericVector stprob;
    vector<double > stprob;
    int seq=0;
// #pragma omp parallel private(seq)
// {
// #pragma omp for schedule(dynamic)
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
      
      // for(int position = 1;position <= SEQUENCES[seq].Size();++position){
      //   double totalprob=0;
      //   for(int j=0;j<names2indicesinprobarray[i].size();++j){
      //     totalprob += Prob[seq][names2indicesinprobarray[i][j]][position];
      //   }
      //   
      //   stprob.push_back(totalprob);
      // }
      
      
      
      
    }
//}
    out_list.push_back(Rcpp::wrap(stprob),print_names[i]);

  }
  
  return out_list;
  
}



Rcpp::List  Forward_Backward_algorithm::getCoverProbDF(){
  extern parameters PARAMS;
  extern NOMeSeqData SEQUENCES;
  extern DNAbind_obj_vector BINDING_OBJECTS;
  
  
  //Rcpp::List out_list(2 + print_names.size());
  Rcpp::List out_list;
  vector<string > seqnames;
  vector<int > positions;
  // fill seq and pos
  for(int seq=0;seq < SEQUENCES.Size();seq++){
    int firstDatPos = SEQUENCES[seq].firstDatpos;
    int lastDatPos = SEQUENCES[seq].lastDatpos;
//    Rcpp::Rcout<<"FirstPos="<<firstDatPos<<"\tLastPos="<<lastDatPos<<endl;
    for(int position = firstDatPos;position <= lastDatPos;++position){
      seqnames.push_back(SEQUENCES[seq].Name());
      positions.push_back(position - firstDatPos + 1);
    }
    //just for testing. uncomment the above code later
    // for(int position = 1;position <= SEQUENCES[seq].Size();++position){
    //   seqnames.push_back(SEQUENCES[seq].Name());
    //   positions.push_back(position - firstDatPos + 1);
    // }
    
  }
  
  //out_list[0] = Rcpp::wrap(seqnames);
  //out_list[1] = Rcpp::wrap(positions);
  out_list.push_back(Rcpp::wrap(seqnames),"seq");
  out_list.push_back(Rcpp::wrap(positions),"pos");
  //out_list.push_back(Rcpp::as<Rcpp::IntegerVector>(positions),"pos");
  
  // fill coverage probabilities

  for(int i=0;i<print_names.size();++i){
    //Rcpp::NumericVector stprob;
    vector<double > covprob;
//#pragma omp for schedule(dynamic)
    for(int seq=0;seq < SEQUENCES.Size();seq++){
      int firstDatPos = SEQUENCES[seq].firstDatpos;
      int lastDatPos = SEQUENCES[seq].lastDatpos;
      
      for(int position = firstDatPos;position <= lastDatPos;++position){
        double coverprob=get_coverage_prob_at_pos_name_index(seq,position,i);
        covprob.push_back(coverprob);
      }
      // just for testing. uncomment the above
      // for(int position = 1;position <= SEQUENCES[seq].Size();++position){
      //   double coverprob=get_coverage_prob_at_pos_name_index(seq,position,i);
      //   covprob.push_back(coverprob);
      // }
      
    }
    out_list.push_back(Rcpp::wrap(covprob),print_names[i]);
    
  }
  
  return out_list;
  
}


Rcpp::List  Forward_Backward_algorithm::getGenomeSummaryDF(){
  
  // calculate statistics across the whole amplicon (all reads)
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
  
  
  // extern parameters PARAMS;
  // string filename = PARAMS.priorfile;
  // if(!filename.empty()){
  //   SetGenomeSummary();
  //   ofstream outputstream(filename.c_str());
  //   outputstream<<"#Name\tPrior\tExp_Numb_Sites\tCoverage\tSite<0.5\tSite>0.5\tPositions<0.5\tPositions>0.5\n";
  //   for(int wm=0;wm<print_names.size();++wm){
  //     outputstream<<print_names[wm]<<"\t"<<genomesummary[wm][0]<<"\t"<<genomesummary[wm][1]<<"\t"<<genomesummary[wm][2]<<"\t"<<genomesummary[wm][3]<<"\t"<<genomesummary[wm][4]<<"\t"<<genomesummary[wm][5]<<"\t"<<genomesummary[wm][6]<<endl;
  //   }
  //   outputstream.close();
  // }
  return out_list;
}



void Forward_Backward_algorithm::clear(){

  F.clear();//swap(vector<vector<double> >());
  R.clear();//(vector<vector<double> >());
  Prob.clear();//(vector<vector<vector<double> > >());
  genomesummary.clear();//(vector<vector<double > > ()); //this vector contains expected prior, expected number of sites and expected coverage for the whole genome and for each TF
  //Forward_Backward_algorithm(){};
  print_indexes.clear();//swap(vector<int >());	// this array contains indexes in object vector that will be printed, i.e. map i - index in Prob array to j - index in object array
  names2indexes.clear();//swap(vector<vector<int > >()); // this array contains map: i - index in print_names to subarray of indexes in object vector with this name (given that for the same tf we create two object with + and - orientation)
  print_names.clear();//swap(vector<string >()); // this array contain names of the objects that will be printed
  
  names2indicesinprobarray.clear();//swap(vector<vector<int > >()); // this array contains map i - index in names to subarray of indices in Prob array
  
}



/*
void Forward_Backward_algorithm::PrintProfilesToFiles()
{
  extern parameters PARAMS;
  extern NOMeSeqData SEQUENCES;

  int first_dat_pos = SEQUENCES.firstDatPos;
  if(!PARAMS.profilefile.empty()){
   Rcpp::Rcout<<PARAMS.profilefile<<endl;
   if(PARAMS.profilefile == "STDOUT"){
    Rcpp::Rcout<<"#=== Begining of the profile of start probabilities ===\n";
    Rcpp::Rcout<<"#Sequence\tPosition";
    for(int i=0;i<print_names.size();++i){
	Rcpp::Rcout <<"\t"<<print_names[i];
    }
    Rcpp::Rcout<<endl;
    for(int seq=0;seq < SEQUENCES.Size();seq++){
	for(int position = 1;position <= SEQUENCES[seq].Size();++position){
		Rcpp::Rcout<<SEQUENCES[seq].Name()<<"\t"<<position - first_dat_pos;
		for(int i=0;i<print_names.size();++i){
			double totalprob=0;
			for(int j=0;j<names2indicesinprobarray[i].size();++j){
				totalprob += Prob[seq][names2indicesinprobarray[i][j]][position];
			}
			Rcpp::Rcout<<"\t"<<scientific<<totalprob;
		}
		Rcpp::Rcout<<endl;
	}
    }
   }
   else{
    ofstream outputstream(PARAMS.profilefile.c_str());
    outputstream<<"#Sequence\tPosition";
    for(int i=0;i<print_names.size();++i){
	outputstream <<"\t"<<print_names[i];
    }
    outputstream<<endl;
    for(int seq=0;seq < SEQUENCES.Size();seq++){
	for(int position = 1;position <= SEQUENCES[seq].Size();++position){
		outputstream<<SEQUENCES[seq].Name()<<"\t"<<position- first_dat_pos;
		for(int i=0;i<print_names.size();++i){
			double totalprob=0;
			for(int j=0;j<names2indicesinprobarray[i].size();++j){
				totalprob += Prob[seq][names2indicesinprobarray[i][j]][position];
			}
			outputstream<<"\t"<<scientific<<totalprob;
		}
		outputstream<<endl;
	}	
    }
    outputstream.close();
   }
  }
}
*/

/*
void Forward_Backward_algorithm::PrintCoverageProfilesToFiles()
{
  
  extern parameters PARAMS;
  extern NOMeSeqData SEQUENCES;
	int first_dat_pos = SEQUENCES.firstDatPos;
  if(!PARAMS.coverageprofile.empty()){
	ofstream outputstream(PARAMS.coverageprofile.c_str());
    	outputstream<<"#Sequence\tPosition";
	for(int i=0;i<print_names.size();++i){
		
		outputstream<<"\t"<<print_names[i];
		
	}
	outputstream<<endl;
	for(int seq=0;seq< SEQUENCES.Size();seq++){
		for(int position = 1;position <= SEQUENCES[seq].Size();++position){
			outputstream<<SEQUENCES[seq].Name()<<"\t"<<position- first_dat_pos;
			for(int i=0;i<print_names.size();++i){
				outputstream<<"\t"<<scientific<<get_coverage_prob_at_pos_name_index(seq,position,i);
			}
			outputstream<<endl;
		}
    	}
	outputstream.close();
  }
  
}

void Forward_Backward_algorithm::PrintGenomeSummaryToFile()
{
	extern parameters PARAMS;
	string filename = PARAMS.priorfile;
	if(!filename.empty()){
		SetGenomeSummary();
		ofstream outputstream(filename.c_str());
		outputstream<<"#Name\tPrior\tExp_Numb_Sites\tCoverage\tSite<0.5\tSite>0.5\tPositions<0.5\tPositions>0.5\n";
		for(int wm=0;wm<print_names.size();++wm){
			outputstream<<print_names[wm]<<"\t"<<genomesummary[wm][0]<<"\t"<<genomesummary[wm][1]<<"\t"<<genomesummary[wm][2]<<"\t"<<genomesummary[wm][3]<<"\t"<<genomesummary[wm][4]<<"\t"<<genomesummary[wm][5]<<"\t"<<genomesummary[wm][6]<<endl;
		}
		outputstream.close();
	}
}

*/

