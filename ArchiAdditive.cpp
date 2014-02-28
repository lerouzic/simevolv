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



#include "ArchiAdditive.h"
#include "Parconst.h"
#include "main.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;



// constructors and destuctor

ArchiAdditive::ArchiAdditive()
{
    assert(false); // The default constructor should never be used.
}


ArchiAdditive::ArchiAdditive(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


ArchiAdditive::ArchiAdditive(const ParameterSet& param)
    : Architecture(param)
{
	// Nothing to do here.
	// mutrate and mutsd are already intialized in the constructor of the parent class
}


// operator overload
/*
ostream& operator << (ostream& out, const ArchiAdditive& archi)
{
    out << "=== Type of model ===" << endl;
    out << endl << "Addititive model" << endl;
    out << endl;
    out << endl;

    out << "=== Genetic map ===" << endl;
    out << archi.gmap;
    out << endl;

    out << "=== Mutation rates ===" << endl;
    out << endl;
    for (int i = 0; i < archi.nb_loc(); i++)
    {
        out << "Loc" << i+1 << "\t" << archi.mutation_rate(i) << endl;
    }
    out << endl;
    return(out);
}
*/


// functions

Phenotype ArchiAdditive::phenotypic_value (const Genotype& genotype) const
{
    unsigned int nloc = nb_loc();
    unsigned int sall = all_size();
    vector<double> sumall_father(nloc);
    vector<double> sumall_mother(nloc);
    vector<double> sumloc(nloc);
    double phenotype=0.0;

    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {
        for (unsigned int all = 0 ; all < sall ; all++)
        {
            sumall_father[loc] += genotype.gam_father.haplotype[loc]->allele[all];
            sumall_mother[loc] += genotype.gam_mother.haplotype[loc]->allele[all];

            sumloc[loc] = sumall_father[loc] + sumall_mother[loc];
        }
    }

    for (unsigned int loc = 0; loc < nloc; loc++)
    {
        phenotype += sumloc[loc];
    }

    return(Phenotype(phenotype));
}

