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
    std::vector<rate_type> mutrates    (archi->mutation_rates(*this));
    std::vector<rate_type> mutmutrates (archi->mutmutation_rates());
    
    for (unsigned int loc = 0; loc < archi->nb_loc(); loc++)
    {
        if (Random::randnum() < mutrates[loc])
        // This is not really efficient (many drawings with low probabilities)
        {
			make_mutation(loc);     
        }
        if ((mutmutrates[loc] > 0.0) && (Random::randnum() < mutmutrates[loc])) {
			make_mutmutation(loc);
		}
    }
}


/* force to make a mutation at a random locus 
 * then > make_mutation(loc) */
void Haplotype::make_mutation(bool test /* = false */)
{ // We have to normalize based on the locus mutation rates
	std::vector<rate_type> mutrates (Architecture::Get()->mutation_rates(*this));
	rate_type sum_mutrates = 0.0;
	for (auto mu  : mutrates) sum_mutrates += mu;
	
	// weighted sampling of the locus
	auto rr = Random::randnum() / sum_mutrates;
	for (size_t loc = 0; loc < nb_loc(); loc++) {
		if (rr < mutrates[loc]) {
			make_mutation(loc, test);
			return;
		}
	}
	assert (false && "Thou shalt not be here!");
}

/* force to make a mutation at a chosen locus :
 * replace the locus (vector of allele) by a new one (<allele_mutation)*/
void Haplotype::make_mutation(unsigned int loc, bool test /* = false */)
{
	Architecture * archi = Architecture::Get();
	shared_ptr<Allele> a = archi->allele_mutation(haplotype[loc], loc, test);
	haplotype[loc] = a;
}

void Haplotype::make_mutmutation()
{
	int loc = floor(Random::randnum()*nb_loc());
	make_mutmutation(loc);
}

void Haplotype::make_mutmutation(unsigned int loc)
{
	Architecture * archi = Architecture::Get();	
	shared_ptr<Allele> a = archi->allele_mut_mutation(haplotype[loc], loc);
	haplotype[loc] = a;
}

std::shared_ptr<const Allele> Haplotype::allele_at_loc(unsigned int loc) const
{
	return(haplotype[loc]);
}

/* perform the recombination between the the father and the mother haplotype */
Haplotype Haplotype::recombine(const Haplotype & hap1, const Haplotype & hap2) 
{
    Architecture * archi = Architecture::Get();

    Haplotype result(hap1);
    bool copy1 = Random::randnum() < 0.5 ? true : false;
    unsigned int nloc = hap1.nb_loc();

    for (unsigned int locus = 0; locus < nloc; locus++)
    {
		if (!copy1) 
		{
			result.haplotype[locus] = hap2.haplotype[locus];
		}
		if ((locus < nloc - 1) && (Random::randnum() < archi -> recombination_rate(locus)))
        { 
            copy1 = !copy1;
        }
	}

    return(result);
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
