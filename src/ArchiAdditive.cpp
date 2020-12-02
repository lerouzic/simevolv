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
#include "Environment.h"

#include "Parconst.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(ArchiAdditive)
#endif

using namespace std;

// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiAdditive::ArchiAdditive(const ParameterSet& param)
    : Architecture(param)
{
    update_param_internal(param); 
}

ArchiAdditive::~ArchiAdditive()
{
}


// inherited functions

/* calculate the phenotypic function depending on the genotype 
 * here : sum of the genotypic values */
Phenotype ArchiAdditive::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
{
	// Epigenetic transmission not implemented yet!
    PhenoVector phenotype(sall);

	for (unsigned int all = 0; all < sall; all++)
		phenotype[all] = 0.0;

    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {
		PhenoVector sumloc = genotype.combine_at_loc(loc, &Allele::combine_add);
		
		for (unsigned int all = 0; all < sall; all++)
			phenotype[all] += sumloc[all];
    }
    
    for (unsigned int all = 0; all < sall; all++)
		if (envir)
			phenotype[all] += Environment::final_disturb();

	Phenotype ans(phenotype);
	ans.scale_transform(transfo);
    return ans;
}

