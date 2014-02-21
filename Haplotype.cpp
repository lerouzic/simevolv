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



#include "Haplotype.h"
#include "Architecture.h"
#include "Parameters.h"
#include "Random.h"

#include <numeric>
#include <string>
#include <cmath>

using namespace std;



// constructors/destructor

Haplotype::Haplotype()
{
    int nloc = Haplotype::nb_loc();
    
    cerr << "Calling Haplotype default constructor: Should probably not happen."<< endl;
    
    for(int i = 0; i < nloc; i++)
    {
		shared_ptr<Allele> a (new Allele());
        haplotype.push_back(a);
    }
}

Haplotype::Haplotype(const ParameterSet & param)
{
    int nloc = Haplotype::nb_loc();
    
    for(int i = 0; i < nloc; i++)
    {
		shared_ptr<Allele> a (new Allele(param));
        haplotype.push_back(a);
    }
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

int Haplotype::nb_loc() const
{
    Architecture * archi = Architecture::Get();
    int nloc = archi -> nb_loc();

    return nloc;
}


void Haplotype::draw_mutation()
{
    Architecture * archi = Architecture::Get();
    for (int loc = 0; loc < nb_loc(); loc++)
    {
        if (Random::randnum() < archi->mutation_rate(loc))
        // This is not really efficient (many drawings with low probabilities)
        {
			make_mutation(loc);        
        }
    }
}


void Haplotype::make_mutation()
{
    int loc = floor(Random::randnum()*nb_loc()); // static_cast<double>(nb_loc())
	make_mutation(loc);
}

void Haplotype::make_mutation(unsigned int loc)
{
	shared_ptr<Allele> a(new Allele(*haplotype[loc]));
	a->make_mutation(loc);
    haplotype[loc] = a;  
}


// output and debug

void Haplotype::write_debug(ostream & out) const
{
    //~ for (vector<Allele>::const_iterator it = haplotype.begin(); it != haplotype.end(); it++)
    //~ {
//~ //        out << *it << "\t";
    //~ }
    out << endl;
}


void Haplotype::write_xml(ostream & out) const
{
    out << "xml output: not implemented yet.\n";
}


void Haplotype::write_simple(ostream& out) const
{
    out << "simple output: not implemented yet.\n";
}


