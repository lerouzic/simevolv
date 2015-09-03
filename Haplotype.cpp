// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
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



#include "Haplotype.h"

#include "Architecture.h"
#include "Parameters.h"
#include "Random.h"

#include <numeric>
#include <string>
#include <cmath>

using namespace std;



// constructors/destructor

/* constructor using the parameters given by the parameters file and the Architecture files*/
Haplotype::Haplotype(const ParameterSet & param)
{
    int nloc = Haplotype::nb_loc();
    Architecture * archi = Architecture::Get();
    
    for(int loc = 0; loc < nloc; loc++)
    {
		shared_ptr<Allele> a = archi->allele_init(param, loc);
        haplotype.push_back(a);
    }
}

/* constructor using a vector of Alleles */
Haplotype::Haplotype(const vector<shared_ptr<Allele>> & vectalleles)
	: haplotype(vectalleles)
{	
}

Haplotype::Haplotype(const Haplotype & templ)
	: haplotype(templ.haplotype)
{	
}


// operator overload

int Haplotype::operator==(const Haplotype& other) const
{
    return((*this).haplotype == other.haplotype);
}

int Haplotype::operator!=(const Haplotype& other) const
{
    return(!(*this == other));
}


// functions

/* collect the number of loci from the Architecture and return it */
unsigned int Haplotype::nb_loc() const
{
    Architecture * archi = Architecture::Get();
    return(archi -> nb_loc());
}

/* determines if there will be a mutation, for each locus of the haplotype 
 * (depending on the mutation rate) 
 * then > make_mutation(loc) */
void Haplotype::draw_mutation()
{
    Architecture * archi = Architecture::Get();
    for (unsigned int loc = 0; loc < archi->nb_loc(); loc++)
    {
        if (Random::randnum() < archi->mutation_rate(loc))
        // This is not really efficient (many drawings with low probabilities)
        {
			make_mutation(loc);     
        }
    }
}


/* force to make a mutation at a random locus 
 * then > make_mutation(loc) */
void Haplotype::make_mutation(bool test /* = false */)
{
    int loc = floor(Random::randnum()*nb_loc()); // static_cast<double>(nb_loc())
	make_mutation(loc, test);
}

/* force to make a mutation at a chosen locus :
 * replace the locus (vector of allele) by a new one (<allele_mutation)*/
void Haplotype::make_mutation(unsigned int loc, bool test /* = false */)
{
	Architecture * archi = Architecture::Get();
	if (test) {
		shared_ptr<Allele> a = archi->allele_mutation_test(haplotype[loc], loc);
		haplotype[loc] = a;		
	} else {
		shared_ptr<Allele> a = archi->allele_mutation(haplotype[loc], loc);
		haplotype[loc] = a;
	}
}

string Haplotype::write_debug() const
{
	ostringstream o;
	for (unsigned int i = 0; i < haplotype.size(); i++) {
		o << "All" << i << ": qq";
		for (unsigned int j = 0; j < haplotype[i]->allele.size(); j++) {
			o << haplotype[i]->allele[j] << " ";
		}
	}
	return(o.str());
}
