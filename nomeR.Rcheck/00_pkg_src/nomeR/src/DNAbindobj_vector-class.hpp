#ifndef _dnabindobjvec_h_
#define _dnabindobjvec_h_
#include "DNAbinding_object-class.hpp"
#include <Rcpp.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include "parameters-class.hpp"
#include <fstream>
#include "binding_object_model-class.hpp"
#include "background-class.hpp"
#include "utils_globvars.hpp"

using namespace std;


class DNAbind_obj_vector{
	vector<DNAbinding_object *> objvector;
	int size; // number of footprint models
public:

	int maxwmlen; // length of longest footprint model
	vector<string > groups; // vector of unique groups present in the input footprint models
	unordered_map<string, vector<int >> group2indices; // map: group -> vector of indices in objvector


	size_t Size() const;
	DNAbind_obj_vector();
	~DNAbind_obj_vector();

	DNAbind_obj_vector(const Rcpp::List _bind_objs,
                    const parameters &params);
	int create(const Rcpp::List _bind_objs,
            const parameters &params);

	const DNAbinding_object* operator [](size_t i) const;

	size_t getGroupsSize() const{
		return groups.size();
	};

	const vector<int >& getGroupIndexVec(string groupName) const{
		return group2indices.at(groupName);
	}

	void print();

	void clear();

	vector<vector<double > > calc_theor_joint_prob(vector<double > ftp_cover_priors, // here vector of priors also represent lengths, namely ith element of the vector
                                                // has length i+1, e.g. ftp_cover_priors[0] is a prior for bg with length 1
                                                double bg_protect_prob,
                                                double footprint_protect_prob,
                                                int max_spacing);

	// method to calculate scores for all footprints, including background given a sequence
	vector<vector<double >> getFtpModelScores(const fragProtectData& fragData) const;


};

#endif
