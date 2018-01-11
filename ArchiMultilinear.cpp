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



#include "ArchiMultilinear.h"
#include "Environment.h"

#include "Parconst.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

BOOST_CLASS_EXPORT(ArchiMultilinear)
#endif

using namespace std;


// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiMultilinear::ArchiMultilinear(const ParameterSet& param)
    : Architecture(param)
    , epsilon2(vector<vector<allele_type>>(0))
    , epsilon3(vector<vector<vector<allele_type>>>(0))
{
	update_param_internal(param); 

    
    // building epsilon matrices
    for (unsigned int loc1 = 0; loc1 < nloc; loc1++)
    {
        if (nloc > 1)
        {
            for (unsigned int loc2 = loc1+1; loc2 < nloc; loc2++)
            {
                set_epsilon2(loc1, loc2, param.getpar(GENET_EPSILON2)->GetDouble());
                if (nloc > 2)
                {
                    for (unsigned int loc3 = loc2+1; loc3 < nloc; loc3++)
                    {
                        set_epsilon3(loc1, loc2, loc3, param.getpar(GENET_EPSILON3)->GetDouble());
                    }
                }
            }
        }
    }
    flag_epistasis2 = !param.getpar(GENET_EPSILON2)->is_nil();
    flag_epistasis3 = !param.getpar(GENET_EPSILON3)->is_nil();
}

ArchiMultilinear::~ArchiMultilinear()
{
}

// functions

/* return the value for 2nd-order epistasis */
allele_type ArchiMultilinear::get_epsilon2(unsigned int loc1, unsigned int loc2) const
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert (loc1 != loc2);
    assert (loc2 < nloc);
    // assert(int(epsilon2.size()) >= loc1+1);
    // assert(int(epsilon2[loc1].size()) >= loc2-loc1);
    return(epsilon2[loc1][loc2-loc1-1]);
}

/* return the value for 3rd-order epistasis */
allele_type ArchiMultilinear::get_epsilon3(unsigned int loc1, unsigned int loc2, unsigned int loc3) const
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    if (loc2 > loc3)
    {
        swap(loc2, loc3);
    }
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert(loc1 != loc2);
    assert(loc2 != loc3);
    assert(loc3 < nloc);
    //assert(int(epsilon3.size()) >= loc1+1);
    //assert(int(epsilon3[loc1].size()) >= loc2-loc1);
    //assert(int(epsilon3[loc1][loc2-loc1-1].size()) >= loc3-loc2);
    return(epsilon3[loc1][loc2-loc1-1][loc3-loc2-1]);
}

/* sets the value of 2nd-order epistasis */
void ArchiMultilinear::set_epsilon2(unsigned int loc1, unsigned int loc2, allele_type value)
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert (loc1 != loc2);
    assert (loc2 < nb_loc());

    while(epsilon2.size() < (loc1+1))
    {
        vector<allele_type> v(0);
        epsilon2.push_back(v);
    }

    while(epsilon2[loc1].size() < (loc2-loc1))
    {
        epsilon2[loc1].push_back(0);
    }

    epsilon2[loc1][loc2-loc1-1] = value;
}

/* sets the value of 3rd-order epistasis */
void ArchiMultilinear::set_epsilon3(unsigned int loc1, unsigned int loc2, unsigned int loc3, allele_type value)
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    if (loc2 > loc3)
    {
        swap(loc2, loc3);
    }
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert(loc1 != loc2);
    assert(loc2 != loc3);
    assert(loc3 < nb_loc());

    while(epsilon3.size() < nb_loc()-2)
    {
        vector<vector<allele_type> > v(0);
        epsilon3.push_back(v);
    }

    while(epsilon3[loc1].size() < (nb_loc()-loc1-2))
    {
        vector<allele_type> v(0);
        epsilon3[loc1].push_back(v);
    }

    while(epsilon3[loc1][loc2-loc1-1].size() < (nb_loc()-loc2-1))
    {
        epsilon3[loc1][loc2-loc1-1].push_back(0);
    }

    epsilon3[loc1][loc2-loc1-1][loc3-loc2-1] = value;
}

/* calculate the phenotypic function depending on the genotype 
	here : sum of the genotypic values, correlation with epistasis values */
Phenotype ArchiMultilinear::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
{
	// Epigenetics not implemented!
    vector<pheno_type> sumloc(nloc);
    pheno_type phenotype = 0.0;

    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {	
		vector<pheno_type> tmp_sum = genotype.combine_at_loc(loc, &Allele::combine_add);

        sumloc[loc] = 0.0;
        for (unsigned int all = 0; all < sall; all++)
			sumloc[loc] += tmp_sum[all];
    }

    for (unsigned int loc1 = 0; loc1 < nloc; loc1++)
    {
        phenotype += sumloc[loc1];
        for (unsigned int loc2 = loc1+1; loc2 < nloc; loc2++)
        {
            if(is_epistasis2())
            {
                phenotype += get_epsilon2(loc1, loc2) * sumloc[loc1] * sumloc[loc2];
            }

            for (unsigned int loc3 = loc2+1; loc3 < nloc; loc3++)
            {
                if(is_epistasis3())
                {
                    phenotype += get_epsilon3(loc1, loc2, loc3) * sumloc[loc1] * sumloc[loc2] * sumloc[loc3];
                }
            }
        }
    }
    if (envir) 
		phenotype += Environment::final_disturb();
    return Phenotype(PhenoScalar(phenotype));
}

