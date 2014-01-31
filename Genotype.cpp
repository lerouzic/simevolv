// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

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

Genotype::Genotype()
    : gam_father(Haplotype())
    , gam_mother(Haplotype())
{
}


Genotype::Genotype(const Haplotype& father, const Haplotype& mother)
    : gam_father(father)
    , gam_mother(mother)
{
}


Genotype::Genotype(const Genotype& copy)
    : gam_father(copy.gam_father)
    , gam_mother(copy.gam_mother)
{
}

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

int Genotype::nb_loc() const
{
    Architecture * archi = Architecture::Get();
    int nloc = archi -> nb_loc();

    return nloc;
}


int Genotype::all_size() const
{
    Architecture * archi = Architecture::Get();
    int sall = archi -> all_size();

    return sall;
}


Haplotype Genotype::recombine() const
{
    Architecture * archi = Architecture::Get();

    const Haplotype * current = NULL;
    const Haplotype * other = NULL;
    Haplotype result;

    int nloc = nb_loc();

    if (Random::randnum() < 0.5)
    {
        current = &gam_father;
        other   = &gam_mother;
    }
    else
    {
        current = &gam_mother;
        other   = &gam_father;
    }

    result.haplotype[0] = current->haplotype[0];

    for (int locus = 1; locus < nloc; locus++)
    {
        if (Random::randnum() < archi -> recombination_rate(locus-1))
        {
            swap(current, other);
        }
        result.haplotype[locus] = current->haplotype[locus];
    }

    return(result);
}


void Genotype::draw_mutation()
{
    gam_father.draw_mutation();
    gam_mother.draw_mutation();
}


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





