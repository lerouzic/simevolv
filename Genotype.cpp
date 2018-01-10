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



#include "Genotype.h"

#include "Random.h"
#include "Architecture.h"
#include "Parameters.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Genotype)
BOOST_CLASS_EXPORT(DiploGenotype)
BOOST_CLASS_EXPORT(HaploGenotype)
#endif

using namespace std;

// functions

/* collect the number of loci from the Architecture and return it */
unsigned int Genotype::nb_loc() const
{
    return (Architecture::Get() -> nb_loc());
}

/* collect the size of the allele from the Architecture and return it */
unsigned int Genotype::all_size() const
{
    return (Architecture::Get() -> all_size());
}

/******************* DiploGenotype *****************************/

/* constructor using two haplotypes */
DiploGenotype::DiploGenotype(const Haplotype& father, const Haplotype& mother)
    : gam_father(father)
    , gam_mother(mother)
{
}

/* copy constructor */
DiploGenotype::DiploGenotype(const DiploGenotype& copy)
    : gam_father(copy.gam_father)
    , gam_mother(copy.gam_mother)
{
}
  
/* constructor using the parameters for building the two haplotypes */
DiploGenotype::DiploGenotype(const ParameterSet & param)
	: gam_father(param)
	, gam_mother(param)
{
}

DiploGenotype* DiploGenotype::clone() const 
{
	return (new DiploGenotype(*this));
}

/* determine if a mutation will occur in one of the haplotype 
 * (depending on the mutation rate) */
void DiploGenotype::draw_mutation()
{
    gam_father.draw_mutation();
    gam_mother.draw_mutation();
}


/* force to make a mutation in one of the haplotype */
void DiploGenotype::make_mutation(bool test /*=false*/)
{
    if (Random::randnum() < 0.5)
    {
        gam_father.make_mutation(test);
    }
    else
    {
        gam_mother.make_mutation(test);
    }
}

Haplotype DiploGenotype::produce_gamete() const 
{
	return(Haplotype::recombine(gam_father, gam_mother));
}

std::vector<allele_type> DiploGenotype::combine_at_loc(unsigned int loc, std::vector<allele_type> (*combineFUN)(const Allele &, const Allele &)) const 
{
	assert(loc < nb_loc());
	return(combineFUN(*gam_father.haplotype[loc], *gam_mother.haplotype[loc]));
}

/******************* HaploGenotype *****************************/

/* constructor using two haplotypes */
HaploGenotype::HaploGenotype(const Haplotype& father, const Haplotype& mother)
    : gam(Haplotype::recombine(father, mother))
{
}

/* copy constructor */
HaploGenotype::HaploGenotype(const HaploGenotype& copy)
    : gam(copy.gam)
{
}
  
/* constructor using the parameters for building the two haplotypes */
HaploGenotype::HaploGenotype(const ParameterSet & param)
	: gam(param)
{
}

HaploGenotype* HaploGenotype::clone() const
{
	return (new HaploGenotype(*this));
}

/* determine if a mutation will occur in one of the haplotype 
 * (depending on the mutation rate) */
void HaploGenotype::draw_mutation()
{
    gam.draw_mutation();
}


/* force to make a mutation in one of the haplotype */
void HaploGenotype::make_mutation(bool test /*=false*/)
{
    gam.make_mutation(test);
}

Haplotype HaploGenotype::produce_gamete() const 
{
	return(gam);
}

std::vector<allele_type> HaploGenotype::combine_at_loc(unsigned int loc, std::vector<allele_type> (*combineFUN)(const Allele &, const Allele &)) const 
{
	assert(loc < nb_loc());
	// For haploid genomes, the "combine" function does not make sense. 
	return(gam.haplotype[loc]->get_raw());
}
