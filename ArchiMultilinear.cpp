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



void Epsilon2::set(size_t i, size_t j, allele_type v)
{
	if (i > j)
		swap(i, j);
	assert(i != j);
	if (fixed_dim) {
		assert(j < dim);
	} else {
		if (j >= dim) dim = j+1;
	}
	
	size_t D1index = i * (dim - 1) + j; // suboptimal
	if (value.size() <= D1index)
		value.resize(D1index+1, 0.0);
	value[D1index] = v;
}

allele_type Epsilon2::get(size_t i, size_t j) const
{
	if (i > j)
		swap(i, j);
	assert (j < dim);
	size_t D1index = i * (dim - 1) + j;
	return(value.at(D1index));
}


void Epsilon3::set(size_t i, size_t j, size_t k, allele_type v)
{
	if (i > j)
		swap(i, j);
	if (i > k)
		swap(i, k);
	if (j > k)
		swap(j, k);

	assert(i < j);
	assert(j < k);
	
	if (fixed_dim) {
		assert(k < dim);
	} else {
		if (k >= dim) dim = k+1;
	}
	
	size_t D1index = i * (dim - 1) * (dim - 1) + j * (dim - 1) + k; // suboptimal
	if (value.size() <= D1index)
		value.resize(D1index+1, 0.0);
	value[D1index] = v;
}

allele_type Epsilon3::get(size_t i, size_t j, size_t k) const
{
	if (i > j)
		swap(i, j);
	if (i > k)
		swap(i, k);
	if (j > k)
		swap(j, k);

	assert(i < j);
	assert(j < k);
	assert(k < dim);
	
	size_t D1index = i * (dim - 1) * (dim - 1) + j * (dim - 1) + k;
	return(value.at(D1index));
}



// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiMultilinear::ArchiMultilinear(const ParameterSet& param)
    : Architecture(param)
    , nphen(param.getpar(GENET_NBPHEN) -> GetInt())    
    , epsilon2(vector<Epsilon2>())
//    , epsilon3(vector<Epsilon3>())
{
	sall = nb_phen();
	update_param_internal(param); 

	assert(nphen > 0);

	for (unsigned int a = 0; a < nphen; a++) {
		for (unsigned int b = 0; b < nphen; b++) {
			for (unsigned int c = 0; c < nphen; c++) {
				// 2-order epistasis
				Epsilon2 eps(nloc);
				if (nloc > 1) { // otherwise, no epistasis possible
					for (unsigned int i = 0; i < nloc-1; i++) {
						for (unsigned int j = i+1; j < nloc; j++) {
							if ((a == b) && (b == c)) {
								eps.set(i, j, param.getpar(GENET_EPSILON2e)->GetDouble());
							} else {
								eps.set(i, j, param.getpar(GENET_EPSILON2p)->GetDouble());
							}
						}
					}
				}
				epsilon2.push_back(eps);
			}
		}
	}
}

// methods


/* calculate the phenotypic function depending on the genotype 
	here : sum of the genotypic values, correlation with epistasis values */
Phenotype ArchiMultilinear::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
{
	// Epigenetics not implemented!
	
	// Step 1: build a vector of genotypes (for each locus)
	vector<PhenoVector> sumv(nloc);
	for (size_t loc = 0; loc < nloc; loc++)
		sumv[loc] = genotype.combine_at_loc(loc, &Allele::combine_add);
	// sumv[locus][trait]
	
	// Step2: compute the phenotype
	// This algorithm features all second-order terms for the multivariate characters multilinear model
	// from Hansen & Wagner 2001. 
	PhenoVector phen(nphen);
	
	for (size_t p = 0; p < nphen; p++) {
		pheno_type tmpsum = 0.0;
		
		// additive part
		for (size_t loc = 0; loc < nloc; loc++)
			tmpsum += sumv[loc][p];
			
		// bilinear epistasis
		for (size_t loc1 = 0; loc1 < nloc-1; loc1++)
		{
			if (nloc > 1)
			for (size_t loc2 = loc1+1; loc2 < nloc; loc2++)
			{
				for (size_t trait1 = 0; trait1 < nphen; trait1++)
				{
					for (size_t trait2 = 0; trait2 < nphen; trait2++)
					{ // there is probably a way to change the loop order in order to call get_epsilon only once per trait combination
						const Epsilon2 & eps = get_epsilon2(p, trait1, trait2); 
						tmpsum += sumv[loc1][trait1] * sumv[loc2][trait2] * eps.get(loc1, loc2);
					}
				}
			}
		}
		phen[p] = tmpsum;
	}
	return Phenotype(phen);
}

const Epsilon2 & ArchiMultilinear::get_epsilon2 (size_t t1, size_t t2, size_t t3) const
{
	size_t index = nphen*nphen*t1 + nphen*t2 + t3;
	return(epsilon2.at(index));
}
