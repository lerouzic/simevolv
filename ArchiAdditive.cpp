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



#include "ArchiAdditive.h"

#include "Parconst.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;



// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiAdditive::ArchiAdditive(const ParameterSet& param)
    : Architecture(param)
{
	// Nothing to do here.
	// mutrate and mutsd are already intialized in the constructor of the parent class
}


// inherited functions

/* calculate the phenotypic function depending on the genotype 
 * here : sum of the genotypic values */
Phenotype ArchiAdditive::phenotypic_value (const Genotype& genotype) const
{
    vector<double> phenotype(sall);

	for (unsigned int all = 0; all < sall; all++)
		phenotype[all] = 0.0;

    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {
		vector<double> sumloc = Allele::combine_add
			(*genotype.gam_father.haplotype[loc], *genotype.gam_mother.haplotype[loc]);
		for (unsigned int all = 0; all < sall; all++)
			phenotype[all] += sumloc[all];
    }

    return(Phenotype(phenotype));
}
