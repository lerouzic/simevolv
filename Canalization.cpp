// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Canalization.h"

#include "Statistics.h"
#include "Parconst.h"

#include "Individual.h"
#include "Population.h"

#include <cassert>
#include <vector>
#include <iomanip>

using namespace std;

Canalization::Canalization(unsigned int nb_tests, const Population & pop)
{
	phen_ready = false;
	fit_ready = false;
	if (nb_tests > 0) {
		for (unsigned int i = 0; i < pop.size(); i++) {
			const Individual & ref = pop.pop[i];
			reference_indiv(ref);
			for (unsigned int test = 0; test < nb_tests; test++) {
				mutant_indiv(ref.test_canalization(1, pop)); // So far: only one mutation per mutant
			}
		} 
	}
}


void Canalization::reference_indiv(Individual ind)
{
	reference.push_back(ind);
	vector<Individual> tmp;
	mutants.push_back(tmp);
}

void Canalization::mutant_indiv(Individual ind) 
{
	assert(!(phen_ready || fit_ready));
	assert(!mutants.empty());
	mutants[mutants.size()-1].push_back(ind);
}

void Canalization::process() 
{
	assert(!mutants.empty());
	assert(!reference.empty());
	
	if (!phen_ready) {
		process_phen();
	}	
	if (!fit_ready) {
		process_fit();
	}	
}

void Canalization::process_phen()
{
	assert(!mutants.empty());
	assert(!reference.empty());	// probably unnecessary
	
	for (unsigned int i = 0; i < mutants.size(); i++) { // individual # i
		vector<Phenotype> data_i;
		for (unsigned int j = 0; j < mutants[i].size(); j++) { // mutant # j
			data_i.push_back(mutants[i][j].get_genot_value());
		}
		PhenotypeStat stat_i(data_i);
		
		mean_per_indiv.push_back(stat_i.means_phen());
		var_per_indiv.push_back(stat_i.vars_phen());
	}	
	
	PhenotypeStat stat_m(mean_per_indiv);
	PhenotypeStat stat_v(var_per_indiv);
	
	mean_of_mean = stat_m.means_phen();
	var_of_mean  = stat_m.vars_phen();
	mean_of_var  = stat_v.means_phen();
	var_of_var   = stat_v.vars_phen();
	
	phen_ready = true;
}

void Canalization::process_fit() 
{
	assert(!mutants.empty());
	assert(!reference.empty());	
	
	vector<vector<double> > dat;
	
	for (unsigned int i = 0; i < mutants.size(); i++) {
		vector<double> tmp;
		for (unsigned int j = 0; j < mutants[i].size(); j++) {
			tmp.push_back(mutants[i][j].get_fitness());
		}
		dat.push_back(tmp);
	}
	MultivariateStat stat_fit(dat);
	
	indiv_fitness_mean = stat_fit.means();
	indiv_fitness_var = stat_fit.vars();
	
	fit_ready = true;
}

Phenotype Canalization::phen_canalization()
{
	if (!phen_ready) {
		process_phen();
	}
	
	return(mean_of_var);
}

double Canalization::fitness_canalization()
{
	if (!fit_ready) {
		process_fit();
	}
	
	UnivariateStat st(indiv_fitness_var);
	return(st.mean());
}
