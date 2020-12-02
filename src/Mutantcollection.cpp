// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Mutantcollection.h"

#include "Fitness.h"

#include <iostream>
#include <vector>

using namespace std;



Mutantcollection::Mutantcollection(unsigned int tests, const Individual & ref, const Population & pop)
	: reference(ref)
{
	Fitness::update(pop); // not sure this is strictly necessary...
	
	for (unsigned int i = 0; i < tests; i++) 
	{
		MiniIndividual idv = ref.test_canalization(1, pop);
		collection.push_back(idv);
	}
};

Mutantcollection::~Mutantcollection()
{ 
}

Mutantcollection & Mutantcollection::operator=(const Mutantcollection & tmpl) 
{
	if (this != &tmpl) 
	{
		reference=tmpl.reference;
		collection = tmpl.collection;
	}
	return(*this);
}

Phenotype Mutantcollection::mean_phen() const
{
    vector<Phenotype> pheno_tmp;
	for (const auto & i : collection) {
		pheno_tmp.push_back(i.phen);
	}
	return Phenotype::mean(pheno_tmp);
}

Phenotype Mutantcollection::var_phen() const
{
    vector<Phenotype> pheno_tmp;
	for (const auto & i : collection) {
		pheno_tmp.push_back(i.phen);
	}
	return Phenotype::var(pheno_tmp);
}

fitness_type Mutantcollection::mean_fit() const
{
	vector<fitness_type> vecd;
	for (const auto & i : collection) 
		vecd.push_back(i.fitness);
	auto fitstat = UnivariateStat<fitness_type>(vecd);
    return fitstat.mean();
}

fitness_type Mutantcollection::var_fit() const
{
	vector<fitness_type> vecd;
	for (const auto & i : collection) 
		vecd.push_back(i.fitness);
	auto fitstat = UnivariateStat<fitness_type>(vecd);
    return fitstat.var();
}

//////////////////// DoubleMutantcollection

DoubleMutantcollection::DoubleMutantcollection(unsigned int tests1, unsigned int tests2, const Individual & reference, const Population & pop) 
{
	Fitness::update(pop); // not sure this is strictly necessary...
	
	for (unsigned int i = 0; i < tests1; i++) 
	{
		const Individual mutant = reference.test_canalization(1, pop);
		Mutantcollection mc(tests2, mutant, pop);
		dcollection.push_back(mc);
	}
}

DoubleMutantcollection::~DoubleMutantcollection() 
{ 
}

DoubleMutantcollection & DoubleMutantcollection::operator=(const DoubleMutantcollection & tmpl)
{
	if (this != &tmpl) 
	{
		dcollection = tmpl.dcollection;
	}
	return(*this);
}

Phenotype DoubleMutantcollection::ref_mean_phen() const 
{
    return Phenotype::mean(ref_phen());
}

Phenotype DoubleMutantcollection::ref_var_phen() const
{
    return Phenotype::var(ref_phen());
}

fitness_type DoubleMutantcollection::ref_mean_fit() const
{
    auto fitstat = UnivariateStat<fitness_type>(ref_fit());
    return fitstat.mean();
}

fitness_type DoubleMutantcollection::ref_var_fit() const
{
    auto fitstat = UnivariateStat<fitness_type>(ref_fit());
    return fitstat.var();
}

vector<Phenotype> DoubleMutantcollection::ref_phen() const
{
	vector<Phenotype> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) 
	{
		ans.push_back(dcollection[i].reference.phen); 
	}
	return(ans);
}

vector<Phenotype> DoubleMutantcollection::var_mutant_phen() const
{
	vector<Phenotype> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) 
	{
		ans.push_back(dcollection[i].var_phen());
	}
	return(ans);
}

vector<fitness_type> DoubleMutantcollection::ref_fit() const
{
	vector<fitness_type> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) 
	{
		ans.push_back(dcollection[i].reference.fitness);
	}
	return(ans);
}

vector<fitness_type> DoubleMutantcollection::var_mutant_fit() const
{
	vector<fitness_type> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) 
	{
		ans.push_back(dcollection[i].var_fit());
	}
	return(ans);
}

