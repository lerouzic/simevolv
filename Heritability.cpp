// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Heritability.h"

#include "Phenotype.h"
#include "Population.h"
#include "Individual.h"
#include "Statistics.h"

#include <vector>

using namespace std;



// constructor

/* Most of the computationally-demanding operations are run in the constructor 
  (simulation of parent-offspring pairs) */
Heritability::Heritability(unsigned int nb_pairs, const Population & pop) 
{
	// nb_pairs is the number of pairs per individual
	
	vector<fitness_type> cumul_equal_fit; 

	// This is a trick to reuse the existing algorithm for sampling individuals uniformly: just assume 
	// a constant fitness function, and generate the corresponding cumulative fitness vector. 
	fitness_type equal_fit = 1.0/pop.size();
	for (unsigned int i = 0; i < pop.size(); i++) 
	{
		cumul_equal_fit.push_back(equal_fit);
		equal_fit = equal_fit + 1.0/pop.size();
	}
	
	for (unsigned int pp = 0; pp < nb_pairs*pop.size(); pp++) 
	{
		/* The algorithm to compute heritability is not very complex. It basically 
		   follows the pattern of a parent - offspring study. 
		   * Step 1: sample two parents randomly in the population
		   * Step 2: cross both parents to get an offspring
		   (here, a database is filled, containing parental and offspring phenotypes and fitnesses)
		   * Step 3 (not in this function): compute parent-offspring regressions and heritabilities. */
		   
		const Individual & Father = pop.pick_parent(cumul_equal_fit);
		const Individual & Mother = pop.pick_parent(cumul_equal_fit);
		Individual Offspring = Individual::mate(Father, Mother);
		Offspring.update_fitness(pop);
		
		ParentOffspring result;
		
		result.father_phen = Father.get_phenotype();
		result.mother_phen = Mother.get_phenotype();
		result.offspring_phen = Offspring.get_phenotype();
		result.father_gen = Father.get_phenotype();
		result.mother_gen = Mother.get_phenotype();
		result.offspring_gen = Offspring.get_phenotype();
		result.father_fit = Father.get_fitness();
		result.mother_fit = Mother.get_fitness();
		result.offspring_fit = Offspring.get_fitness();
		
		parentoffspring.push_back(result);
	}
}


// functions

Phenotype Heritability::h2() const 
{
	// This function calculates the narrow-sense heritability based on the
	// mid-parent - offspring covariance. 
	// In practice, for each phenotypic dimension, a MultivariateStat object is filled with 
	// mid-parent phenotypic value and offspring phenotype. The slope of the regression between
	// them provides a good estimate of heritability. 
	vector<pheno_type> result;
	for(unsigned int trait = 0; trait < parentoffspring[0].father_phen.dimensionality(); trait++) 
	{
		vector<pheno_type> midpar;
		vector<pheno_type> offspring;
		for (unsigned int i = 0; i < parentoffspring.size(); i++) 
		{
			midpar.push_back(0.5*parentoffspring[i].father_phen[trait] + 0.5*parentoffspring[i].mother_phen[trait]);
			offspring.push_back(parentoffspring[i].offspring_phen[trait]);
		}
		vector<vector<pheno_type> > data;
		data.push_back(midpar);
		data.push_back(offspring);
		MultivariateStat<pheno_type> stat(data);
		result.push_back(stat.regression_slope(1, 0));
	}
	return(Phenotype(result));
}

fitness_type Heritability::fit_h2() const 
{
	// Heritability for fitness. Basically the same algorithm than for phenotypes. 
	// Is it necessary to factorize the code? This could generate more problems than 
	// what it solves. 
	vector<fitness_type> midpar;
	vector<fitness_type> offspring;
	for (unsigned int i = 0; i < parentoffspring.size(); i++) 
	{
		midpar.push_back(0.5*parentoffspring[i].father_fit + 0.5*parentoffspring[i].mother_fit);
			offspring.push_back(parentoffspring[i].offspring_fit);
	}
	vector<vector<fitness_type> > data;
	data.push_back(midpar);
	data.push_back(offspring);
	MultivariateStat<fitness_type> stat(data);
	return(stat.regression_slope(1,0));
}
