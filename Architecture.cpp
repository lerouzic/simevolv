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

Architecture::Architecture()
{
    assert(false); // The default constructor should never be used.
}


Architecture::Architecture(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


Architecture::Architecture(const ParameterSet& param)
    : gmap (param)
    , nloc (param.getpar(GENET_NBLOC) -> GetInt())
    , sall (1)
    , mutrate (vector<double> (0))
    , mutsd (vector<double> (0))
{
    for (int i = 0; i < nloc; i++)
    {
        mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i));
        mutsd.push_back(param.getpar(GENET_MUTSD)->GetDouble(i));
    }

}


// instance and initialization

Architecture* Architecture::instance = NULL;


void Architecture::initialize(const ParameterSet& param)
{
    if (Architecture::instance != NULL)
    {
		cerr << "Strange, the Architecture was already intialized. This should not happen." << endl;
        delete Architecture::instance;
        Architecture::instance = NULL;
    }
    //Architecture::instance = new Architecture(param);

    string type_archi = param.getpar(TYPE_ARCHI)->GetString();
    if (type_archi==AR_add)
    {
        Architecture::instance = new ArchiAdditive(param);
    }
    else if (type_archi==AR_mult)
    {
        Architecture::instance = new ArchiMultilinear(param);
    }
    else if (type_archi==AR_reg)
    {
        Architecture::instance = new ArchiRegulatory(param);
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




// operator overload
/*
ostream& operator << (ostream& out, const Architecture& archi)
{
    out << "=== Genetic map ===" << endl;
    out << archi.gmap;
    out << endl;

    out << "=== Mutation rates ===" << endl;
    for (int i = 0; i < archi.nb_loc(); i++)
    {
        out << "Loc" << i+1 << "\t" << archi.mutation_rate(i) << endl;
    }
    out << endl;
    return(out);
}
*/

// functions

int Architecture::nb_loc() const
{
    return nloc;
}


int Architecture::all_size() const
{
    return sall;
}


double Architecture::mutation_rate(int locus) const
{
    assert (locus >= 0);
    assert (locus < nloc);
    return(mutrate[locus]);
}


double Architecture::mutation_sd(int locus) const
{
    assert(locus >= 0);
    assert(locus < nloc);
    return(mutsd[locus]);
}


double Architecture::recombination_rate(int locus) const
{
    assert(locus >=0);
    assert(locus < nloc-1);
    return(gmap.recombination_rate(locus));
}


// for inheritance (virtual function)

shared_ptr<Allele> Architecture::allele_init(const ParameterSet & param) const 
{
	vector<double> tmp;
	for(int i = 0; i < sall; i++)
    {
        tmp.push_back(param.getpar(INIT_ALLELES) -> GetDouble());
    }
    shared_ptr<Allele> a(new Allele(tmp));
    return(a);
}


shared_ptr<Allele> Architecture::allele_mutation(const Allele & templ, unsigned int loc /* = 0 */) const 
{
	// A mutation affects randomly one of the "sites" of the allele (?)

    int mutated_site = floor(sall*Random::randnum());
    double modifier = mutation_sd(loc) * Random::randgauss();
    shared_ptr<Allele> a(new Allele(templ));
    a->allele[mutated_site] += modifier;
    return(a);
}

