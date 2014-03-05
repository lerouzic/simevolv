// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <vector>

#include "Heritability.h"
#include "Phenotype.h"
#include "Population.h"
#include "Individual.h"

using namespace std;

Heritability::Heritability(unsigned int nb_pairs, const Population & pop) {
	
	vector<double> cumul_equal_fit; 

	// This is a trick to reuse the existing algorithm for sampling individuals uniformly: just assume 
	// a constant fitness function.
	double equal_fit = 1.0/pop.size();
	for (unsigned int i = 0; i < pop.size(); i++) {
		cumul_equal_fit.push_back(equal_fit);
		equal_fit = equal_fit + 1.0/pop.size();
	}
	
	for (unsigned int pp = 0; pp < nb_pairs; pp++) {
		const Individual & Father = pop.pick_parent(cumul_equal_fit);
		const Individual & Mother = pop.pick_parent(cumul_equal_fit);
		Individual Offspring = Individual::mate(Father, Mother);
		Offspring.update_fitness(pop);
		
		ParentOffspring result;
		
		result.father_phen = Father.get_phenotype();
		result.mother_phen = Mother.get_phenotype();
		result.offspring_phen = Offspring.get_phenotype();
		result.father_gen = Father.get_genot_value();
		result.mother_gen = Mother.get_genot_value();
		result.offspring_gen = Offspring.get_genot_value();
		result.father_fit = Father.get_fitness();
		result.mother_fit = Mother.get_fitness();
		result.offspring_fit = Offspring.get_fitness();
		
		parentoffspring.push_back(result);
	}
}

Phenotype Heritability::h2() const {
	vector<double> result;
	for(unsigned int trait = 0; trait < parentoffspring[0].father_phen.dimensionality(); trait++) {
		vector<double> midpar;
		vector<double> offspring;
		for (unsigned int i = 0; i < parentoffspring.size(); i++) {
			midpar.push_back(0.5*parentoffspring[i].father_phen[trait] + 0.5*parentoffspring[i].mother_phen[trait]);
			offspring.push_back(parentoffspring[i].offspring_phen[trait]);
		}
		vector<vector<double> > data;
		data.push_back(midpar);
		data.push_back(offspring);
		MultivariateStat stat(data);
		result.push_back(stat.regression_slope(1, 0));
	}
	return(Phenotype(result));
}

Phenotype Heritability::H2() const {
	// the code is quite duplicated with h2().
	vector<double> result;
	for(unsigned int trait = 0; trait < parentoffspring[0].father_phen.dimensionality(); trait++) {
		vector<double> midpar;
		vector<double> offspring;
		for (unsigned int i = 0; i < parentoffspring.size(); i++) {
			midpar.push_back(0.5*parentoffspring[i].father_gen[trait] + 0.5*parentoffspring[i].mother_gen[trait]);
			offspring.push_back(parentoffspring[i].offspring_gen[trait]);
		}
		vector<vector<double> > data;
		data.push_back(midpar);
		data.push_back(offspring);
		MultivariateStat stat(data);
		result.push_back(stat.regression_slope(1, 0));
	}
	return(Phenotype(result));
}

double Heritability::fit_h2() const {
	vector<double> midpar;
	vector<double> offspring;
	for (unsigned int i = 0; i < parentoffspring.size(); i++) {
		midpar.push_back(0.5*parentoffspring[i].father_fit + 0.5*parentoffspring[i].mother_fit);
			offspring.push_back(parentoffspring[i].offspring_fit);
	}
	vector<vector<double> > data;
	data.push_back(midpar);
	data.push_back(offspring);
	MultivariateStat stat(data);
	return(stat.regression_slope(1,0));
}