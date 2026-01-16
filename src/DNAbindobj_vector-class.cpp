#include "DNAbindobj_vector-class.hpp"

DNAbind_obj_vector::DNAbind_obj_vector(){

}

DNAbind_obj_vector::~DNAbind_obj_vector(){

    for(size_t wm = 0; wm < objvector.size(); wm++){
        if(objvector[wm] != NULL){
            delete objvector[wm];
        }
    }
}

const DNAbinding_object * DNAbind_obj_vector::operator [](size_t i) const{
    if(i>objvector.size()-1 || i<0){
        Rcpp::stop("DNAbind_obj_vector::operator[](size_t i):  The index is out of range");
    }
    return objvector[i];
}

size_t DNAbind_obj_vector::Size() const{
    return size;
}
DNAbind_obj_vector::DNAbind_obj_vector(const Rcpp::List _bind_objs,
                                       const parameters &params){
    create(_bind_objs,
           params);
}

int DNAbind_obj_vector::create(const Rcpp::List _bind_objs,
                               const parameters &params){

    maxwmlen = 1;

    for(size_t wm=0; wm < _bind_objs.size(); ++wm){
        Rcpp::List pinf = Rcpp::as<Rcpp::List >(_bind_objs[wm]);

        if(!pinf.containsElementNamed("PROTECT_PROB")){
            Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element PROTECT_PROB\n";
            return(0);
        }
        if(!pinf.containsElementNamed("PRIOR")){
            Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element PRIOR\n";
            return(0);
        }
        if(!pinf.containsElementNamed("NAME")){
            Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element NAME\n";
            return(0);
        }
        if(!pinf.containsElementNamed("GROUP")){
            Rcpp::Rcerr<<"DNAbind_obj_vector::create: Error! At least one element in list of sequences does not contain element GROUP\n";
            return(0);
        }

        vector<double > prot_prob = Rcpp::as<vector<double > >(pinf["PROTECT_PROB"]);
        double prior = Rcpp::as<double >(pinf["PRIOR"]);
        string name = Rcpp::as<string >(pinf["NAME"]);
        string group = Rcpp::as<string >(pinf["GROUP"]);

        DNAbinding_object *newobj=new binding_object_model(prot_prob,
                                                           prior,
                                                           name,
                                                           group);
        objvector.push_back(newobj);

        // check if group exists in groups and add if not
        if(find(groups.begin(), groups.end(), newobj->group) == groups.end()){
            groups.push_back(newobj->group); // add group into the vector of unique groups
            vector<int > tmp(1,objvector.size() - 1);
            group2indices[newobj->group] = tmp;
        } else{
            group2indices[newobj->group].push_back(objvector.size() - 1);
        }

        if(newobj->len > maxwmlen)
            maxwmlen = newobj->len;

    }


    //initialise backgound model. always at the end of the vector
    Background *bg = new Background(params);
    objvector.push_back(bg);
    // check if background group exists in groups and add if not.
    // sometimes we use short footprints (2-5bp) to account for correlated background noise.
    if(find(groups.begin(), groups.end(), bg->group) == groups.end()){
        groups.push_back(bg->group); // add group into the vector of unique groups
        vector<int > tmp(1,objvector.size() - 1);
        group2indices[bg->group] = tmp;
    } else{
        group2indices[bg->group].push_back(objvector.size() - 1);
    }
    size = objvector.size();


    // normalize priors so that they sum up to 1
    double priorsum=0;
    for(size_t i=0;i<objvector.size();++i){
        priorsum += objvector[i]->prior;// * objvector[i]->len;
    }

    for(int i=0;i<objvector.size();++i){
        objvector[i]->prior = (objvector[i]->prior)/priorsum;
    }

    // calculate posteriors for infinite non-informative sequence

    double sum = 0;
    for(int i = 0; i < objvector.size(); ++i){
        sum += objvector[i]->prior * objvector[i]->len;
    }

    for(int i = 0; i < objvector.size(); ++i){
        objvector[i]->nonInformPosterior = objvector[i]->prior / sum;
    }
    return objvector.size();
}


// method to calculate scores for all footprints, including background given a sequence;
vector<vector<double >> DNAbind_obj_vector::getFtpModelScores(const fragProtectData& fragData) const{
    vector<vector<double >> ftpScoresMatrix;
    for(int wm = 0; wm < size; ++wm){
        ftpScoresMatrix.push_back(objvector[wm]->get_seq_scores_vec(fragData));
    }
    return ftpScoresMatrix;
}



void DNAbind_obj_vector::clear(){

    maxwmlen=0;

    for(size_t wm=0;wm<objvector.size();wm++){
        if(objvector[wm] != NULL){
            delete objvector[wm];
        }
    }
    objvector.clear();
    size=0;
}


void DNAbind_obj_vector::print(){

    for(size_t wm = 0;wm < objvector.size(); ++wm){
        objvector[wm]->print_normalized();
    }

}





// Function to calculate theoretical joint probabiltiies to observe 00,01,10,11 at distance S

vector<vector<double > > DNAbind_obj_vector::calc_theor_joint_prob(vector<double > ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector
                                                                   // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
                                                                   double bg_protect_prob,
                                                                   double footprint_protect_prob,
                                                                   int max_spacing){

    // define vector with ftp lengths
    vector<int > ftp_lengths;
    double R_const=0;
    for(size_t ftp=0;ftp < ftp_cover_priors.size();++ftp){
        ftp_lengths.push_back(ftp + 1);
        R_const += ftp_cover_priors[ftp]/(ftp + 1);
    }

    if(R_const <= 0){
        Rcpp::stop("DNAbind_obj_vector::calc_theor_joint_prob: ERROR! Value of R_const is negative or zero.");
    }
    // calculate start priors for ftps
    vector<double > ftp_start_priors;
    for(size_t ftp=0; ftp < ftp_cover_priors.size();++ftp){
        ftp_start_priors.push_back((ftp_cover_priors[ftp]/ftp_lengths[ftp])/R_const);
    }

    // calculate forward partition sum
    vector<double > FW(max_spacing + 1,0); // here FW[0] is intial condition for FW
    vector<double > ProbBg(max_spacing + 1,0);
    vector<double > sigma_cumul(max_spacing + 1,0);
    FW[0] = 1;
    ProbBg[0] = 0;
    sigma_cumul[0] = 0;
    for(int d=1; d <= max_spacing; ++d){
        FW[d] = 0;
        for(size_t ftp=0; ftp < ftp_start_priors.size(); ++ftp){
            if(d - ftp_lengths[ftp] >=0){
                FW[d] += ftp_start_priors[ftp] * FW[d - ftp_lengths[ftp]];
            }
        }

        // set bg prob
        ProbBg[d] = ftp_start_priors[0] * FW[d-1];
        sigma_cumul[d] = sigma_cumul[d-1] + ProbBg[d];
    }

    // set emission probs vectors for convenience
    double beta1 = bg_protect_prob;
    double beta0 = 1-bg_protect_prob;
    double alpha1 = footprint_protect_prob;
    double alpha0 = 1-footprint_protect_prob;
    vector<double > alpha1_vec;
    vector<double > alpha0_vec;
    alpha1_vec.push_back(beta1);
    alpha0_vec.push_back(beta0);
    for(size_t ftp=1; ftp < ftp_start_priors.size(); ++ftp){
        alpha1_vec.push_back(alpha1);
        alpha0_vec.push_back(alpha0);
    }

    // calculate joint theoretical probabilities
    vector<vector<double > > Pjoint;
    // fill at S=1, i.e probs of 1 and 0. no distance
    vector<double > tmp(5,0);
    tmp[0] = 1;
    tmp[1] = alpha0 + (beta0 - alpha0) * ftp_cover_priors[0];
    tmp[2] = 0;
    tmp[3] = 0;
    tmp[4] = alpha1 + (beta1 - alpha1) * ftp_cover_priors[0];
    Pjoint.push_back(tmp);
    for(int S = 2; S <= max_spacing; ++S){
        vector<double > prob_joint(5,0); // This correspond to Spacing, P(0,0|S), P(0,1|S), P(1,0|S), P(1,1|S)
        double sigma_sum0 = 0;
        double sigma_sum1 = 0;
        for(int ftp=0; ftp < ftp_start_priors.size(); ++ftp){
            int ftp_len = ftp_lengths[ftp];
            if(S - ftp_len - 1 >0){
                sigma_sum0 += ftp_start_priors[ftp] * alpha0_vec[ftp] * sigma_cumul[S - ftp_len - 1];
                sigma_sum1 += ftp_start_priors[ftp] * alpha1_vec[ftp] * sigma_cumul[S - ftp_len - 1];
            }
        }

        prob_joint[0] = S;
        // calc probs

        // P(0,0|S)
        prob_joint[1] = alpha0 * alpha0 + R_const * ftp_start_priors[0] * alpha0 * (beta0 - alpha0) +
            R_const * (beta0 - alpha0) * ( (alpha0 + (beta0 - alpha0) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum0);
        // P(0,1|S)
        prob_joint[2] = alpha1 * alpha0 + R_const * ftp_start_priors[0] * alpha1 * (beta0 - alpha0) +
            R_const * (beta1 - alpha1) * ( (alpha0 + (beta0 - alpha0) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum0);
        // P(1,0|S)
        prob_joint[3] = alpha0 * alpha1 + R_const * ftp_start_priors[0] * alpha0 * (beta1 - alpha1) +
            R_const * (beta0 - alpha0) * ( (alpha1 + (beta1 - alpha1) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum1);
        // P(1,1|S)
        prob_joint[4] = alpha1 * alpha1 + R_const * ftp_start_priors[0] * alpha1 * (beta1 - alpha1) +
            R_const * (beta1 - alpha1) * ( (alpha1 + (beta1 - alpha1) * ftp_start_priors[0]) * sigma_cumul[S-1] - sigma_sum1);


        Pjoint.push_back(prob_joint);
    }

    return(Pjoint);

}

