// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
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
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;



// constructors/destructor

/* default constructor  -  should never be used */
Architecture::Architecture()
{
    assert(false);
}


/* copy constructor  -  should never be used */
Architecture::Architecture(const Architecture& archi)
{
    assert(false); 
}

/* constructor using the paramater given by ParameterSet */
Architecture::Architecture(const ParameterSet& param)
    : gmap (param)
    , nloc (param.getpar(GENET_NBLOC) -> GetInt())
    , sall (1) // by default, change for the architecture that needed a vector as allele
    , mutrate (vector<double> (0))
    , mutsd (vector<double> (0))
{
    for (unsigned int i = 0; i < nloc; i++)
    {
        mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i));
        mutsd.push_back(param.getpar(GENET_MUTSD)->GetDouble(i));
    }

}


// instance and initialization

/* put the existence of the architecture to non-existent */
Architecture* Architecture::instance = NULL;


/* initialization of the global architecture of the genetic system,
 * depending on the architecture type */
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
    else if (type_archi==AR_masel)
    {
        Architecture::instance = new ArchiMasel(param);
    }
    else if (type_archi==AR_siegal)
    {
        Architecture::instance = new ArchiSiegal(param);
    }
    else
    {
        assert("Wrong architecture type");
    }
}


Architecture* Architecture::Get()
{
    assert(Architecture::instance != NULL && "Attempt to access Architecture::Get() before initialization.");
    return(Architecture::instance);
}


// functions

/* return the number of loci */
unsigned int Architecture::nb_loc() const
{
    return nloc;
}


/* return the size of the allele */
unsigned int Architecture::all_size() const
{
    return sall;
}


/* return the mutation rate at a given locus */
double Architecture::mutation_rate(unsigned int locus) const
{
    assert (locus >= 0);
    assert (locus < nloc);
    return(mutrate[locus]);
}


/* return the mutation effect at a given locus */
double Architecture::mutation_sd(unsigned int locus) const
{
    assert(locus >= 0);
    assert(locus < nloc);
    return(mutsd[locus]);
}


/* return the recombination rate at a given locus */
double Architecture::recombination_rate(unsigned int locus) const
{
    assert(locus >=0);
    assert(locus < nloc-1);
    return(gmap.recombination_rate(locus));
}


// for inheritance (virtual function)

/* initialization of the alleles (1 allele = 1 vector of value) */
shared_ptr<Allele> Architecture::allele_init(const ParameterSet & param, unsigned int loc /* = 0 */) const 
{
	// Here we don't need to know the locus, but inherited classes may.
	vector<double> tmp;
	for(unsigned int i = 0; i < sall; i++)
    {
        tmp.push_back(param.getpar(INIT_ALLELES) -> GetDouble());
    }
    shared_ptr<Allele> a(new Allele(tmp));
    return(a);
}


/* force to make a mutation at one position of the allele :
 * replace the value at the mutated site by a new value */
shared_ptr<Allele> Architecture::allele_mutation(const Allele & templ, unsigned int loc /* = 0 */) const 
{
    int mutated_site = floor(sall*Random::randnum());
    double modifier = mutation_sd(loc) * Random::randgauss();
    shared_ptr<Allele> a(new Allele(templ));
    a->allele[mutated_site] += modifier;
    return(a);
}

