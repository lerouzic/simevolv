// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>
// Copyright 2015      Christine Mayer     <christine.mayer@ibv.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "ArchiBoolean.h"

#include "Parconst.h"
#include "Random.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

BOOST_CLASS_EXPORT(ArchiBoolean)
#endif

using namespace std;


// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiBoolean::ArchiBoolean(const ParameterSet& param)
: Architecture(param)
, bucket_matrix(vector<vector<float_type>>(0))
, logic_operator(vector<float_type>(0))
, ploc (param.getpar(PHEN_NBLOC) -> GetInt())
, type(param.getpar(SCALE)->GetString())

{
    update_param_internal(param); 

    
    float_type threshold_matrix = param.getpar(MATRIX_DENS)->GetDouble();
    float_type threshold_operator = param.getpar(LOG_OPERATOR_DENS)->GetDouble();
    
    // building matrices
    //TODO: Test for "empty" lines -> if a line has no ones, then the locus in the phenotype is not build by any locus of the genotype!
    vector<float_type> v(0);
    for (unsigned int loc1 = 0; loc1 < ploc; loc1++)
    {
        for (unsigned int loc2 = 0; loc2 < nloc; loc2++ )
        {
            if( Random::randnum()< threshold_matrix)
                v.push_back(1);
            else
                v.push_back(0);
        }
        set_bucket_matrix(loc1, v);
        v.erase(v.begin(),v.end());
        
    }
    
    //set logic_operator vector - this vector codes for the logical operation that are used per phenotype locus
    for(unsigned int loc1 = 0; loc1< ploc; loc1++)
    {
        if(Random::randnum() < threshold_operator)
            logic_operator.push_back(1);
        else
            logic_operator.push_back(0);
    }
}

ArchiBoolean::~ArchiBoolean()
{
}

// functions

unsigned int ArchiBoolean::nb_phen() const {
    return ploc; // Not sure (ALR 5 Juil 2017)
}

/* return value of bucket-matrix */
float_type ArchiBoolean::get_bucket_matrix(unsigned int loc1, unsigned int loc2) const
{
    return(bucket_matrix[loc1][loc2]);
}

//return value of the vector logic_operator
float_type ArchiBoolean::get_logic_operator(unsigned int loc) const
{
    return(logic_operator[loc]);
}

/* sets row of bucket_matrix */
//TODO: should set whole matrix --> get the setting out of the constructor
void ArchiBoolean::set_bucket_matrix(unsigned int loc1, vector<float_type> value)
{
    bucket_matrix.push_back(value);
}

//TODO: setter for logical_operator, to get it out of constructor
//TODO: Implementation of print Matrix
//TODO: print logical_operator
//TODO: print genotype
//TODO: print phenotype
//string ArchiBoolean::print_epsilon2() const
//{
//    
//}

//init of the haploid genotype -> gets random number of gaussian distribution. Parmeters (mean, sd) are defined in the parameter file.
//TODO: Maybe change to uniform distribution? (Random:randnum) -> then parameter has to be changed! Parameter_gaussian to Parameter_double
shared_ptr<Allele> ArchiBoolean::allele_init(const ParameterSet & param, unsigned int loc /* = 0 */) const
{
    vector<float_type> tmp;
    
    if(param.getpar(INIT_ALLELES) -> GetDouble()<0.5)
        tmp.push_back(0);
    else
        tmp.push_back(1);
    
    
    shared_ptr<Allele> a;
    a = shared_ptr<Allele>(new Allele(tmp));
   
    return(a);
}


/* calculate the phenotypic function depending on the genotype
	here : boolean logic operations between loci of the genotype */
Phenotype ArchiBoolean::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
	// Warning: the variable "envir" has no effect in this model. 
	// Warning: epigenetic transmission not implemented yet.
{
    pheno_type phenotype_calc =1;
    vector<pheno_type> phenotype(0);
    int phenotypesum =0;
    
//calculation of the diploid genotype, based on OR operation. Each locus of the father and the mother gets combined by an OR function
//function is defined in Allele.cpp
    vector<vector<allele_type>> y;
    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {
		y.push_back(genotype.combine_at_loc(loc, &Allele::combine_OR));
    }

    
//Implementation of AND and OR using the multilinear model (formulas)
//AND: z=y1*y2*....*yn
//OR: (By definition it is a NOT AND) z=1-(1-y1)(1-y2)...(1-yn)
    for(unsigned int locp = 0;locp<ploc; locp++){
        phenotype_calc=1;
        // 0 == AND || 1 == OR
        if(logic_operator[locp]==0){
            for(unsigned int locg = 0; locg<nloc;locg++){
                if(bucket_matrix[locp][locg]==1){
                    phenotype_calc = phenotype_calc*y[locg][0];
                }
            }
            phenotype.push_back(phenotype_calc);
        }
        else if(logic_operator[locp]==1){
            for(unsigned int locg = 0; locg<nloc;locg++){
                if(bucket_matrix[locp][locg]==1){
                    phenotype_calc = phenotype_calc*(1-y[locg][0]);
                }
            }
            phenotype.push_back(1-phenotype_calc);
        }
        
    }
    
//return of the phenotype or calcualtion of phenotype. Depending on the parameter SCALE in the parameter file.
// 1. SC_vector returns the whole phenotype vector
// 2. SC_int returns the sum of the phenotype vector.
// 3. SC_dec treats the vector as binary number and converts it to decimal
// 4. SC_combi is a combination of 1. and 2. The first number in the vector is the sum of the vector, followed by the whole vector
    if(type == SC_vector){
        return Phenotype(PhenoVector(phenotype));
    }
    else if(type==SC_int){
    
        for(unsigned int i=0;i<phenotype.size();i++)
        {
            phenotypesum=phenotypesum + phenotype[i];
        }
        //cout<< phenotypesum << endl;
        return Phenotype(PhenoScalar(phenotypesum));
    }
    else if(type==SC_dec){
        for(unsigned int i=0;i<phenotype.size();i++)
        {
            phenotypesum=phenotypesum + (phenotype[phenotype.size()-i-1]*pow(2,i));
        }
        //cout<< phenotypesum << endl;
        return Phenotype(PhenoScalar(phenotypesum));
    }
    else if(type==SC_combi){
      
        for(unsigned int i=0;i<phenotype.size();i++)
        {
            phenotypesum=phenotypesum + phenotype[i];
        }
        phenotype.insert(phenotype.begin(),phenotypesum);
        
        return Phenotype(PhenoVector(phenotype));
    }

    assert(type==SC_vector||type==SC_dec||type==SC_combi||type==SC_int);
    return Phenotype(); // should never happen
}
