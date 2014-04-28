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



#include "Genotype.h"
#include "Random.h"
#include "Architecture.h"
#include "Parameters.h"

#include <algorithm>
#include <cassert>
#include <cmath>


using namespace std;


// constructors and destuctor

/* default constructor */
Genotype::Genotype()
    : gam_father(Haplotype())
    , gam_mother(Haplotype())
{
}


/* constructor using two haplotypes */
Genotype::Genotype(const Haplotype& father, const Haplotype& mother)
    : gam_father(father)
    , gam_mother(mother)
{
}


/* copy constructor */
Genotype::Genotype(const Genotype& copy)
    : gam_father(copy.gam_father)
    , gam_mother(copy.gam_mother)
{
}
 
 
/* constructor using the parameters for building the two haplotypes */
Genotype::Genotype(const ParameterSet & param)
	: gam_father(param)
	, gam_mother(param)
{
}


// operator overload

int Genotype::operator== (const Genotype& other) const
{
    return
    (
        ((this->gam_father == other.gam_father) &&
        (this->gam_mother == other.gam_mother))
        ||
        ((this->gam_father == other.gam_mother) &&
        (this->gam_mother == other.gam_father))
    );
}


int Genotype::operator!= (const Genotype& other) const
{
    return(!(*this==other));
}


// functions

/* collect the number of loci from the Architecture and return it */
int Genotype::nb_loc() const
{
    Architecture * archi = Architecture::Get();
    int nloc = archi -> nb_loc();
    return nloc;
}


/* collect the size of the allele from the Architecture and return it */
int Genotype::all_size() const
{
    Architecture * archi = Architecture::Get();
    int sall = archi -> all_size();

    return sall;
}


/* perform the recombination between the the father and the mother haplotype */
Haplotype Genotype::recombine() const
{
    Architecture * archi = Architecture::Get();

    Haplotype result(gam_father);
    bool copy_father = Random::randnum() < 0.5 ? true : false;
    unsigned int nloc = nb_loc();

    for (unsigned int locus = 0; locus < nloc; locus++)
    {
		if (!copy_father) {
			result.haplotype[locus] = gam_mother.haplotype[locus];
		}
		if ((locus < nloc - 1) && (Random::randnum() < archi -> recombination_rate(locus)))
        { 
            copy_father = !copy_father;
        }
	}

    return(result);
}


/* determine if a mutation will occur in one of the haplotype 
 * (depending on the mutation rate) */
void Genotype::draw_mutation()
{
    gam_father.draw_mutation();
    gam_mother.draw_mutation();
}


/* force to make a mutation in one of the haplotype */
void Genotype::make_mutation()
{
    int gen = floor(2*Random::randnum()+1);
    if (gen==1)
    {
        gam_father.make_mutation();
    }
    else
    {
        gam_mother.make_mutation();
    }
}


// output

void Genotype::write_debug(ostream & out) const
{
    out << "Gamete 1" << endl;
    gam_father.write_debug(out);
    out << "Gamete 2" << endl;
    gam_mother.write_debug(out);
}


void Genotype::write_xml(ostream & out) const
{
    out << "xml output: not implemented yet.\n";
}


void Genotype::write_simple(ostream& out) const
{
    out << "simple output: not implemented yet.\n";
}





