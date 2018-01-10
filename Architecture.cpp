// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Architecture.h"

#include "ArchiAdditive.h"
#include "ArchiMultilinear.h"
#include "ArchiRegulatoryMatrix.h"
#include "ArchiBoolean.h"
#include "Parconst.h"
#include "Random.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(Architecture)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(Architecture)
#endif

using namespace std;

// constructors/destructor

/* constructor using the paramater given by ParameterSet */
Architecture::Architecture(const ParameterSet& param)
    : gmap (param)                                    // This cannot change when the parameter file changes
    , nloc (param.getpar(GENET_NBLOC) -> GetInt())    // The same, this cannot change even if the parameter file changes
    , sall (1) // by default, change for the architecture that needed a vector as allele
    , mutrate (vector<rate_type> (0))
    , mutsd (vector<allele_type> (0))
    , mutsd_test (vector<allele_type> (0))
    , plasticity_strength (vector<pheno_type> (0))
    , plasticity_signal (vector<pheno_type> (0))
    , iofile ("")
{
	// update_param_internal(param); // This should be called in derived classes only
       
    // This happens only when the param constructor is called.
    if (param.exists(FILE_ARCHI))
		iofile = param.getpar(FILE_ARCHI) -> GetString();
}

Architecture::~Architecture() 
{
}

// instance and initialization

/* put the existence of the architecture to non-existent */
Architecture* Architecture::instance = NULL;

/* initialization of the global architecture of the genetic system, depending on the architecture type */
void Architecture::initialize(const ParameterSet& param)
{
    if (Architecture::instance != NULL)
    {
		cerr << "Strange, the Architecture was already intialized. This should not happen." << endl;
        delete Architecture::instance;
        Architecture::instance = NULL;
    }

	string type_archi = param.getpar(TYPE_ARCHI)->GetString();

    if (type_archi==AR_add)
    {
		Architecture::instance = new ArchiAdditive(param);
    }
    else if (type_archi==AR_mult)
    {		
		Architecture::instance = new ArchiMultilinear(param);
    }
    else if (type_archi==AR_wagner)
    {
        Architecture::instance = new ArchiWagner(param);
    }
    else if (type_archi==AR_siegal)
    {
        Architecture::instance = new ArchiSiegal(param);
    }
    else if (type_archi==AR_m2)
    {
		Architecture::instance = new ArchiM2(param);
    }
    else if(type_archi==AR_Boolean)
        Architecture::instance = new ArchiBoolean(param);
    else
    {
        cerr << "Wrong architecture type" << endl;
        exit(EXIT_FAILURE);
    }

    if (Architecture::instance == NULL) {
		cerr << "Error when creating the genetic architecture." << endl;
		exit(EXIT_FAILURE);
	}
}

/* Reads the serialized architecture file */
void Architecture::initialize (const string & archi_file) 
{	
    if (Architecture::instance != NULL)
    {
		cerr << "Strange, the Architecture was already intialized. This should not happen." << endl;
        delete Architecture::instance;
        Architecture::instance = NULL;
    }
    
    try {
		ifstream infile(archi_file);
		if (infile.good()) {
			boost::archive::text_iarchive ar(infile);
			ar >> Architecture::instance;
		}
	}
	catch (std::exception & e)
	{
		cerr << "Exception " << e.what() << " when reading architecture file " << archi_file << endl;
		exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Error when reading architecture file " << archi_file << endl;
		exit(EXIT_FAILURE);
	}
}

// Updates the parameter values (at least, those that can change meaningfully during simulation)
void Architecture::update_param(const ParameterSet& param)
{
	Architecture * archi = Architecture::Get();
	archi-> update_param_internal(param);
}

Architecture* Architecture::Get()
{
    assert(Architecture::instance != NULL && "Attempt to access Architecture::Get() before initialization.");
    return(Architecture::instance);
}

void Architecture::terminate() 
{
	delete Architecture::instance;
}

// functions

/* return the number of loci */
unsigned int Architecture::nb_loc() const
{
    return nloc;
}

unsigned int Architecture::nb_phen() const
{
    return 1; // not a bad choice for most architectures
}

/* return the size of the allele */
unsigned int Architecture::all_size() const
{
    return sall;
}

/* return the mutation rate at a given locus */
rate_type Architecture::mutation_rate(unsigned int locus) const
{
    assert (locus < nloc);
    return(mutrate[locus]);
}

/* return the mutation effect at a given locus */
allele_type Architecture::mutation_sd(unsigned int locus) const
{
    assert(locus < nloc);
    return(mutsd[locus]);
}

/* return the mutation effect for canalization tests at a given locus */
allele_type Architecture::mutation_sd_test(unsigned int locus) const
{
    assert(locus < nloc);
    return(mutsd_test[locus]);
}

/* return the recombination rate at a given locus */
rate_type Architecture::recombination_rate(unsigned int locus) const
{
    assert(locus < nloc-1);
    return(gmap.recombination_rate(locus));
}


// for inheritance (virtual function)

/* initialization of the alleles (1 allele = 1 vector of value) */
shared_ptr<Allele> Architecture::allele_init(const ParameterSet & param, unsigned int loc /* = 0 */) const 
{
	// Here we don't need to know the locus, but inherited classes may.
	vector<allele_type> tmp;
	for(unsigned int i = 0; i < sall; i++)
    {
        tmp.push_back(param.getpar(INIT_ALLELES) -> GetDouble());
    }
  
    string type_alleles = param.getpar(TYPE_ALLELES) -> GetString();
    shared_ptr<Allele> a;
    if (type_alleles==TA_norm)
    {
        a = shared_ptr<Allele>(new Allele(tmp));
    }
    else if (type_alleles==TA_zero)
    {
		a = shared_ptr<Allele>(new Allele_zero(tmp));
    } 
    else
    {
		cerr << "Unknown allele type -- theoretically this should not happen" << endl;
		exit(EXIT_FAILURE);
	}
    return(a);
}

/* force to make a mutation at one position of the allele : replace the value at the mutated site by a new value */
shared_ptr<Allele> Architecture::allele_mutation(const shared_ptr<Allele> templ, unsigned int loc /* = 0 */) const 
{
    return(templ->make_mutant(mutation_sd(loc)));
}

shared_ptr<Allele> Architecture::allele_mutation_test(const shared_ptr<Allele> templ, unsigned int loc /* = 0 */) const 
{
    return(templ->make_mutant(mutation_sd_test(loc)));
}

/* Updates parameters when the parameter set changes. 
 * Only parameters that are meaningful to change are updated */
void Architecture::update_param_internal(const ParameterSet& param)
{
	mutrate.clear();
	mutsd.clear();
	mutsd_test.clear();
	plasticity_strength.clear();
    plasticity_signal.clear();
    
    for (unsigned int i = 0; i < nloc; i++)
    {
		if (param.getpar(GENET_MUTTYPE)->GetString() == MT_locus)
			mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i));
		else 
			mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i)/static_cast<rate_type>(nloc));
        mutsd.push_back(param.getpar(GENET_MUTSD)->GetDouble(i));
        mutsd_test.push_back(param.getpar(OUT_CANAL_MUTSD)->GetDouble(i));        
    }
    
    if (param.exists(ENVIRO_PLASTICITY)) {
        for (unsigned int i = 0; i < this->nb_phen(); i++) {
            plasticity_strength.push_back(param.getpar(ENVIRO_PLASTICITY)->GetDouble(i));
            plasticity_signal.push_back(param.getpar(FITNESS_OPTIMUM)->GetDouble(i));
        }
    } else {
        plasticity_strength = vector<pheno_type>(this->nb_phen(), 0.0);
        plasticity_signal = vector<pheno_type>(this->nb_phen(), 0.0);
    }

    
}
