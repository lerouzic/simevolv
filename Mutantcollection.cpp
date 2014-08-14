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
	, phenostat(NULL)
	, fitstat(NULL)
{
	Fitness::update(pop); // not sure this is strictly necessary...
	
	for (unsigned int i = 0; i < tests; i++) {
		MiniIndividual idv = ref.test_canalization(1, pop);
		collection.push_back(idv);
	}
};

Mutantcollection::~Mutantcollection()
{ // memory leak if the pointers are not deleted
	delete phenostat;
	delete fitstat;
}

Mutantcollection & Mutantcollection::operator=(const Mutantcollection & tmpl) 
{
	if (this != &tmpl) {
		reference=tmpl.reference;
		collection = tmpl.collection;
		phenostat=NULL;
		fitstat=NULL;
	}
	return(*this);
}

Phenovec Mutantcollection::mean_phen() const
{
	if (phenostat == NULL)
		compute_phenostat();
	return(phenostat->means_phen());
}

Phenovec Mutantcollection::var_phen() const
{
	if (phenostat == NULL)
		compute_phenostat();
	return(phenostat->vars_phen());
}

double Mutantcollection::mean_fit() const
{
	if (fitstat == NULL)
		compute_fitstat();
	return(fitstat->mean());
}

double Mutantcollection::var_fit() const
{
	if (fitstat == NULL)
		compute_fitstat();
	return(fitstat->var());
}

void Mutantcollection::compute_phenostat() const
{ // This function can be const because phenostat is mutable
	if (phenostat != NULL) {
		cerr << "Warning: mutantcollection is updated several times." << endl;
		delete phenostat; // otherwise the memory leaks!
	}
	vector<Phenotype> vecp;
	for (unsigned int i = 0; i < collection.size(); i++) {
		vecp.push_back(collection[i].phen);
	}
	phenostat = new PhenotypeStat(vecp);
}

void Mutantcollection::compute_fitstat() const
{ // this function can be const because fitstat is mutable
	if (fitstat != NULL) {
		cerr << "Warning: mutantcollection is updated several times." << endl;
		delete fitstat; // otherwise the memory leaks!
	}
	vector<double> vecd;
	for (unsigned int i = 0; i < collection.size(); i++) {
		vecd.push_back(collection[i].fitness);
	}
	fitstat = new UnivariateStat(vecd);
}

//////////////////// DoubleMutantcollection

DoubleMutantcollection::DoubleMutantcollection(unsigned int tests1, unsigned int tests2, const Individual & reference, const Population & pop) 
	: refstat(NULL)
	, reffitstat(NULL)
{
	Fitness::update(pop); // not sure this is strictly necessary...
	
	for (unsigned int i = 0; i < tests1; i++) {
		const Individual mutant = reference.test_canalization(1, pop);
		Mutantcollection mc(tests2, mutant, pop);
		dcollection.push_back(mc);
	}
}

DoubleMutantcollection::~DoubleMutantcollection() 
{ 
	delete refstat;
	delete reffitstat;
}

DoubleMutantcollection & DoubleMutantcollection::operator=(const DoubleMutantcollection & tmpl)
{
	if (this != &tmpl) {
		dcollection = tmpl.dcollection;
		// it is probably safer and easier to recompute the statistics stored in cache
		refstat = NULL;
		reffitstat = NULL;
	}
	return(*this);
}

Phenovec DoubleMutantcollection::ref_mean_phen() const 
{
	if (refstat == NULL)
		compute_refstat();
	return(refstat->means_phen());
}

Phenovec DoubleMutantcollection::ref_var_phen() const
{
	if (refstat == NULL)
		compute_refstat();
	return(refstat->vars_phen());
}

double DoubleMutantcollection::ref_mean_fit() const
{
	if (reffitstat == NULL)
		compute_reffitstat();
	return(reffitstat->mean());
}

double DoubleMutantcollection::ref_var_fit() const
{
	if (reffitstat == NULL)
		compute_reffitstat();
	return(reffitstat->var());
}

vector<Phenovec> DoubleMutantcollection::ref_phen() const
{
	vector<Phenovec> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) {
		ans.push_back(dcollection[i].reference.phen.get_pheno()); 
	}
	return(ans);
}

vector<Phenovec> DoubleMutantcollection::var_mutant_phen() const
{
	vector<Phenovec> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) {
		ans.push_back(dcollection[i].var_phen());
	}
	return(ans);
}

vector<double> DoubleMutantcollection::ref_fit() const
{
	vector<double> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) {
		ans.push_back(dcollection[i].reference.fitness);
	}
	return(ans);
}

vector<double> DoubleMutantcollection::var_mutant_fit() const
{
	vector<double> ans;
	for (unsigned int i = 0; i < dcollection.size(); i++) {
		ans.push_back(dcollection[i].var_fit());
	}
	return(ans);
}


void DoubleMutantcollection::compute_refstat() const
{ // this can be const because refstat is mutable
	if (refstat != NULL) {
		cerr << "Warning: DoubleMutantcollection is updated several times." << endl;
		delete refstat; // otherwise the memory leaks!
	}
	refstat = new PhenotypeStat(ref_phen());
}

void DoubleMutantcollection::compute_reffitstat() const
{
	if(reffitstat != NULL) {
		cerr << "Warning: DoubleMutantcollection is updated several times." << endl;
		delete reffitstat; // otherwise the memory leaks!
	}		
	reffitstat = new UnivariateStat(ref_fit());
}
