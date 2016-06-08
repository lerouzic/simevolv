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
#include <cmath> // log function and isnan

using namespace std;

MiniCanIndiv::MiniCanIndiv(const vector<Individual> & variants, const Individual & ref, bool logvar, bool meancentered)
{
	vector<Phenotype> phenos;
	vector<double> fitnesses;
	for (auto variant : variants) {
		phenos.push_back(variant.get_phenotype());
		fitnesses.push_back(variant.get_fitness());
	}
	PhenotypeStat phenostat(phenos);
	UnivariateStat fitnessstat(fitnesses);	
	
	Phenovec rawmeanphen = phenostat.means_phen();
	Phenovec rawvarphen = phenostat.vars_phen();
	double rawmeanfit = fitnessstat.mean_log();
	double rawvarfit = fitnessstat.var_log();
	
	if (!meancentered) {
		Phenovec ref_pheno = ref.get_phenotype().get_pheno();
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
}


/**************************** Canalization ***********************************/

Canalization::Canalization(unsigned int can_tests, const Population & pop, bool logvar, bool meancentered)
{
}

Canalization::~Canalization() { }

Phenovec Canalization::meanpop_canphen() const 
{
	vector<Phenovec> phenpop;
	for (auto minican : popcan)
		phenpop.push_back(minican.canpheno);
	PhenotypeStat phenostat(phenpop);
	return(phenostat.means_phen());
}

Phenovec Canalization::varpop_canphen() const 
{
	vector<Phenovec> phenpop;
	for (auto minican : popcan)
		phenpop.push_back(minican.canpheno);
	PhenotypeStat phenostat(phenpop);
	return(phenostat.vars_phen());
}

double Canalization::meangene_meanpop_canphen() const
{
	vector<double> meanpop;
	for (auto i : meanpop_canphen())
		meanpop.push_back(i);
	UnivariateStat uvstat(meanpop);
	return(uvstat.mean());
}

vector<double> Canalization::meangene_canphen() const
{
	vector<double> ans;
	for (auto minican : popcan) {
		vector<double> indcan;
		for (auto i : minican.canpheno)
			indcan.push_back(i);
		UnivariateStat uvstat(indcan);
		ans.push_back(uvstat.mean());
	}
	return(ans);
}

double Canalization::meanpop_canlogfit() const
{
	UnivariateStat uvstat(canlogfit());
	return(uvstat.mean());
}

double Canalization::varpop_canlogfit() const
{
	UnivariateStat uvstat(canlogfit());
	return(uvstat.var());
}

vector<double> Canalization::canlogfit() const
{
	vector<double> fitnesses;
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
