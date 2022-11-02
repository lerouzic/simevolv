// Copyright 2022       Arnaud Le Rouzic    <arnaud.le-rouzic@universite-paris-saclay.fr>


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "ArchiFKL.h"
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

BOOST_CLASS_EXPORT(ArchiFKL)
#endif

using namespace std;


/* constructor using the paramater given in Architecture and the parameters files */
ArchiFKL::ArchiFKL(const ParameterSet& param)
    : ArchiAdditive(param)
    , nphen(param.getpar(GENET_NBPHEN) -> GetInt())
    , mutmutrate (vector<rate_type> (0))    
{
	assert(nphen > 0);
	sall = nphen;
	update_param_internal(param);
}

unsigned int ArchiFKL::nb_phen() const {
	return(nphen);
}

std::vector<rate_type> ArchiFKL::mutation_rates(const Haplotype & hap) const
{
	std::vector<rate_type> mutrates(nb_loc());
	rate_type tmp_sum = 0.0;
	for (size_t loc = 0; loc < nb_loc(); loc++) {
		const Allele_mut & tmp_all = dynamic_cast<const Allele_mut &>(*hap.allele_at_loc(loc));
		mutrates[loc] = exp(tmp_all.get_mutrate());
		tmp_sum += mutrates[loc];
	}
	for (size_t loc = 0; loc < nb_loc(); loc++) {
		mutrates[loc] *= mutrate[loc]/tmp_sum;
	}
	
	return mutrates;
}

std::vector<rate_type> ArchiFKL::mutmutation_rates() const
{
	return mutmutrate;
}

Phenotype ArchiFKL::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
{
    PhenoVector phenotype(nphen);

	for (unsigned int all = 0; all < nphen; all++)
		phenotype[nphen] = 0.0;

    for (unsigned int loc = 0 ; loc < nloc ; loc++)
    {
		PhenoVector sumloc = genotype.combine_at_loc(loc, &Allele::combine_add);
		assert(sumloc.dimensionality() == nphen + 1); // The first slot of the allele stands for the mutation rate
		
		for (unsigned int all = 0; all < nphen; all++)
			phenotype[all] += sumloc[all+1];
    }
    
    for (unsigned int all = 0; all < nphen; all++)
		if (envir)
			phenotype[all] += Environment::final_disturb();

	Phenotype ans(phenotype);
	ans.scale_transform(transfo);
    return ans;
}

shared_ptr<Allele> ArchiFKL::allele_init(const ParameterSet & param, unsigned int loc /* = 0 */) const 
{  
	vector<string> type_allele_loc;
	vector<allele_type> tmp;
	
	// The first allele slot is the (log) mutation rate
	tmp.push_back(log(mutrate[loc]));
	
	for(unsigned int i = 0; i < sall-1; i++)
    {
		if (param.exists(INIT_ALLELES_FULL))
			tmp.push_back(param.getpar(INIT_ALLELES_FULL) -> GetDouble(loc*nphen+i));
		else
			tmp.push_back(param.getpar(INIT_ALLELES) -> GetDouble());
    }
    
    for (unsigned int i = 0; i < sall; i++)
		type_allele_loc.push_back(param.getpar(TYPE_ALLELES) -> GetString(loc*sall+i));

    shared_ptr<Allele> a(new Allele(tmp, type_allele_loc));

    return(a);
}

/* Replace the value at the mutated site by a new value */
shared_ptr<Allele> ArchiFKL::allele_mutation(const shared_ptr<Allele> templ, unsigned int loc /* = 0 */, bool test /* = false */) const 
{
	if (test) {
		return(templ->make_mutant_all_sites(mutmodels_test[loc]));
	} else {
		return(templ->make_mutant_all_sites(mutmodels[loc]));
	}
}

shared_ptr<Allele_mut> ArchiFKL::allele_mut_mutation(shared_ptr<const Allele_mut> templ, unsigned int loc /* = 0 */) const 
{
	return(templ->make_mutant_mutation(mutmodels[loc]));
}


/* Updates parameters when the parameter set changes. 
 * Only parameters that are meaningful to change are updated */
void ArchiFKL::update_param_internal(const ParameterSet& param)
{
	Architecture::update_param_internal(param);
	
	mutmutrate.clear();
	
	nphen = param.getpar(GENET_NBPHEN) -> GetInt();
	
	for (unsigned int i = 0; i < nloc; i++)
	{
		if (param.getpar(GENET_MUTTYPE)->GetString() == MT_locus)
			mutmutrate.push_back(param.getpar(GENET_MUTMUTRATES)->GetDouble(i));
		else 
			mutmutrate.push_back(param.getpar(GENET_MUTMUTRATES)->GetDouble(i)/static_cast<rate_type>(nloc));
	}
}
