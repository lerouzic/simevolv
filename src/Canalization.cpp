// Copyright 2013-2017      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

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
#include <cmath> // log function and isnan

using namespace std;

MiniCanIndiv::MiniCanIndiv(const vector<Individual> & variants, const Individual & ref, bool logvar, bool meancentered)
{
	vector<Phenotype> phenos;
	vector<fitness_type> fitnesses;
	for (auto variant : variants) {
		phenos.push_back(variant.get_phenotype());
		fitnesses.push_back(variant.get_fitness());
	}

	UnivariateStat<fitness_type> fitnessstat(fitnesses);	
	
	Phenotype rawmeanphen = Phenotype::mean(phenos);
	Phenotype rawvarphen = Phenotype::var(phenos);
	fitness_type rawmeanfit = fitnessstat.mean_log();
	fitness_type rawvarfit = fitnessstat.var_log();
	
	if (!meancentered) {
		Phenotype ref_pheno = ref.get_phenotype();
        // rawvarphen = rawvarphen + (rawmeanphen - ref_pheno)*(rawmeanphen - ref_pheno);
		for (unsigned int i = 0; i < rawmeanphen.dimensionality(); i++) {
			rawvarphen[i] = rawvarphen[i] + (rawmeanphen[i] - ref_pheno[i])*(rawmeanphen[i] - ref_pheno[i]);
		}
		rawvarfit = rawvarfit + (rawmeanfit - ref.get_fitness())*(rawmeanfit - ref.get_fitness());
	}
	
	if (logvar) {
		for (unsigned int i = 0; i < rawmeanphen.dimensionality(); i++) {
			rawvarphen[i] = log(rawvarphen[i]);
			if ((rawvarphen[i] < MIN_LOG_VAR) || (std::isnan(rawvarphen[i])))
				rawvarphen[i] = MIN_LOG_VAR;
		}
		rawvarfit = log(rawvarfit);
		if ((rawvarfit < MIN_LOG_VAR) || (std::isnan(rawvarfit)))
			rawvarfit = MIN_LOG_VAR;
	}
	
	canpheno = rawvarphen;
	canfitness = rawvarfit;
	
	for (unsigned int i = 0; i < rawmeanphen.dimensionality(); i++) {
		vcov.push_back(Phenotype::vcov(phenos, i));
	}
	
}


/**************************** Canalization ***********************************/

Canalization::Canalization(unsigned int can_tests, const Population & pop, bool logvar, bool meancentered)
{
}

Canalization::~Canalization() { }

Phenotype Canalization::meanpop_canphen() const 
{
	vector<Phenotype> phenpop;
	for (auto minican : popcan)
		phenpop.push_back(minican.canpheno);
	return(Phenotype::mean(phenpop));
}

Phenotype Canalization::varpop_canphen() const 
{
	vector<Phenotype> phenpop;
	for (auto minican : popcan)
		phenpop.push_back(minican.canpheno);
	return(Phenotype::var(phenpop));
}

Phenotype Canalization::meanpop_vcov(unsigned int i) const
{
	vector<Phenotype> tmp_vcov_i;
	for (auto minican : popcan)
		tmp_vcov_i.push_back(minican.vcov[i]);
	return(Phenotype::mean(tmp_vcov_i));
}

vector<Phenotype> Canalization::meanpop_vcov() const
{
	vector<Phenotype> mean_vcov;
	const size_t dim_phen = popcan[0].vcov.size(); // Clearly not elegant
	for (unsigned int i = 0; i < dim_phen; i++) {
		mean_vcov.push_back(this->meanpop_vcov(i));
	}
	return(mean_vcov);
}

Phenotype Canalization::varpop_vcov(unsigned int i) const
{
	vector<Phenotype> tmp_vcov_i;
	for (auto minican : popcan)
		tmp_vcov_i.push_back(minican.vcov[i]);
	return(Phenotype::var(tmp_vcov_i));
}

vector<Phenotype> Canalization::varpop_vcov() const
{
	vector<Phenotype> var_vcov;
	const size_t dim_phen = popcan[0].vcov.size(); // Clearly not elegant
	for (unsigned int i = 0; i < dim_phen; i++) {
		var_vcov.push_back(this->varpop_vcov(i));
	}
	return(var_vcov);
}

fitness_type Canalization::meanpop_canlogfit() const
{
	UnivariateStat<fitness_type> uvstat(canlogfit());
	return(uvstat.mean());
}

fitness_type Canalization::varpop_canlogfit() const
{
	UnivariateStat<fitness_type> uvstat(canlogfit());
	return(uvstat.var());
}

vector<fitness_type> Canalization::canlogfit() const
{
	vector<fitness_type> fitnesses;
	for (auto minican : popcan) {
		fitnesses.push_back(minican.canfitness);
	}
	return(fitnesses);
}

/************************ Genetic Canalization ****************************/

GeneticCanalization::GeneticCanalization(unsigned int can_tests, const Population & pop, bool logvar, bool meancentered)
	: Canalization(can_tests, pop, logvar, meancentered)
{
	// In theory, this should not be necessary. In practice, something in the population changes and the Fitness function complains
	Fitness::update(pop);
	
	if (can_tests > 0) 
	{
		for (unsigned int i = 0; i < pop.size(); i++) 
		{
			vector<Individual> mutants;			
			for (unsigned int test = 0; test < can_tests; test++) 
			{
				mutants.push_back(pop.pop[i].test_canalization(1, pop)); 
					// So far: only one mutation per mutant
			}
			MiniCanIndiv minican(mutants, pop.pop[i], logvar, meancentered);
			popcan.push_back(minican);
		}
	}
}

/************************ Init Disturbance Canalization **********************/

DisturbCanalization::DisturbCanalization(unsigned int can_tests, const Population & pop, bool logvar, bool meancentered)
	: Canalization(can_tests, pop, logvar, meancentered)
{
	Fitness::update(pop);
	
	if (can_tests > 0) 
	{
		for (unsigned int i = 0; i < pop.size(); i++) 
		{
			vector<Individual> variants;			
			for (unsigned int test = 0; test < can_tests; test++) 
			{
				variants.push_back(pop.pop[i].test_disturb(pop)); 
			}
			MiniCanIndiv minican(variants, pop.pop[i], logvar, meancentered);
			popcan.push_back(minican);
		}
	}
}

/************************* Environmental Canalization *************************/

EnviroCanalization::EnviroCanalization(unsigned int can_tests, const Population & pop, bool logvar, bool meancentered)
	: Canalization(can_tests, pop, logvar, meancentered)
{
	Fitness::update(pop);
	
	if (can_tests > 0) 
	{
		for (unsigned int i = 0; i < pop.size(); i++) 
		{
			vector<Individual> variants;			
			for (unsigned int test = 0; test < can_tests; test++) 
			{
				variants.push_back(pop.pop[i].test_enviro(pop)); 
			}
			MiniCanIndiv minican(variants, pop.pop[i], logvar, meancentered);
			popcan.push_back(minican);
		}
	}
}
