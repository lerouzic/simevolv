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
#include "Fitness.h"

#include <cassert>
#include <vector>
#include <iomanip>

using namespace std;



// Constructors 

Canalization::Canalization(unsigned int can_tests, const Population & pop)
{
	phen_ready = false;
	fit_ready = false;

	// In theory, this should not be necessary. In practice, something in the population changes and the Fitness function complains
	Fitness::update(pop);
	
	// Fills the database of the object: a collection of reference ("wild") individuals (the individuals of the population), and for each wild individual, nb_tests mutants. 
	if (can_tests > 0) 
	{
		for (unsigned int i = 0; i < pop.size(); i++) 
		{
			const Individual & ref = pop.pop[i];
			reference_indiv(ref);
			for (unsigned int test = 0; test < can_tests; test++) 
			{
				mutant_indiv(ref.test_canalization(1, pop)); // So far: only one mutation per mutant
			}
		} 
	}
}


// User interface.

/* Get the phenotypic canalization scores (i.e. the average across individuals of mutant variances) */
Phenotype Canalization::phen_canalization()
{
	if (!phen_ready) 
	{
		process_phen();
	}
	
	return(mean_of_var);
}

/*	Get the canalization scores for fitness (the average of mutant variances) */
double Canalization::fitness_canalization()
{
	if (!fit_ready) 
	{
		process_fit();
	}
	
	UnivariateStat st(indiv_fitness_var);
	return(st.mean());
}

// Internal functions

/* Sets a new reference individual
Note that, in practice, reference individuals are never used in the calculation. Yet, it is mandatory to provide them, as they indicate that the next mutants will concern another individual.*/ 
void Canalization::reference_indiv(Individual ind)
{
	reference.push_back(ind);
	vector<Individual> tmp;
	mutants.push_back(tmp);
}

/* Adds a new mutant in the database. Note that reference individuals and corresponding mutants have to be entered sequencially, which is not very conveninent (no way to enter first all reference individuals, and then all mutants. Nevermind, this is internal code. */
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
	
	if (!phen_ready) 
	{
		process_phen();
	}	
	
	if (!fit_ready) 
	{
		process_fit();
	}	
}

void Canalization::process_phen()
{
	assert(!mutants.empty());
	assert(!reference.empty());	// probably unnecessary

	// Here is the calculation of the canalization score itself. 
	// The algorithm computes the mean and variance across mutants in each reference individual. 
	// It then computes both the mean and the variance of these means and variances.
	// Everything is stored, and the chosen index is returns by other "getter" functions.
	
	for (unsigned int i = 0; i < mutants.size(); i++) 
	{ // individual # i
		vector<Phenotype> data_i;
		for (unsigned int j = 0; j < mutants[i].size(); j++) 
		{ // mutant # j
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
	// canalization scores for fitness. 
	// The algorithm is very similar to the one for phenotypes, except that fitnesses are unidimensional. 
	assert(!mutants.empty());
	assert(!reference.empty());	
	
	vector<vector<double> > dat;
	
	for (unsigned int i = 0; i < mutants.size(); i++) 
	{
		vector<double> tmp;
		for (unsigned int j = 0; j < mutants[i].size(); j++) 
		{
			tmp.push_back(mutants[i][j].get_fitness());
		}
		dat.push_back(tmp);
	}
	
	MultivariateStat stat_fit(dat);
	
	indiv_fitness_mean = stat_fit.means();
	indiv_fitness_var = stat_fit.vars();
	
	fit_ready = true;
}
