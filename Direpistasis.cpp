// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Direpistasis.h"

#include "Population.h"
#include "Phenotype.h"

#include <vector>
#include <cmath>

using namespace std;



Direpistasis::Direpistasis(unsigned int tests, const Population & pop)
{	
	unsigned int mutants2 = sqrt(tests);
	unsigned int mutants1 = mutants2;
	
	// Below these values, results cannot be calculated. 
	// Just don't call Direpistasis with lower values
	if (mutants1 < 1) mutants1 = 2;
	if (mutants2 < 1) mutants2 = 2;
		
	vector<Phenovec> dir_ind;
	vector<double> fit_dir;
	for (unsigned int i = 0; i < pop.size(); i++) 
	{
		DoubleMutantcollection dmutcol(mutants1, mutants2, pop.pop[i], pop);
		// We can do this because Direpistasis is a friend of Population

		dir_ind.push_back(individual_direpi(dmutcol));
		fit_dir.push_back(individual_fitdir(dmutcol));
	}
	PhenotypeStat dir_stat(dir_ind);
	dir_epi_mean = dir_stat.means_phen();
	dir_epi_var  = dir_stat.vars_phen();
	
	UnivariateStat fit_stat(fit_dir);
	fit_epi_mean = fit_stat.mean();
	fit_epi_var = fit_stat.var();	
}

Phenovec Direpistasis::phen_direpistasis() const
{
	return(dir_epi_mean);
}

Phenovec Direpistasis::var_phen_direpistasis() const
{
	return(dir_epi_var);
}

double Direpistasis::fitness_direpistasis() const
{
	return(fit_epi_mean);
}

double Direpistasis::var_fitness_direpistasis() const
{
	return(fit_epi_var);
}

Phenovec Direpistasis::individual_direpi(const DoubleMutantcollection & dmutcol) const
{ // returns the directional epistasis for a single individual
	
	Phenovec mut1_var = dmutcol.ref_var_phen();
	vector<Phenovec> mut1_phen = dmutcol.ref_phen();
	vector<Phenovec> mut2_var_phen = dmutcol.var_mutant_phen();
	
	vector<double> dir; // directional epistasis
	// loop over characters
	for (unsigned int t = 0; t < mut1_var.dimensionality(); t++) 
	{ 
		
		// Step 1: transpose the results to get means and variances per character
		vector<double> y1_t;
		vector<double> sd_m;
		for (unsigned int i = 0; i < mut1_phen.size(); i++) 
		{
			y1_t.push_back(mut1_phen[i][t]);
			sd_m.push_back(sqrt(mut2_var_phen[i][t]));
		}
		
		// Step 2: compute the regression between the standard deviation of the effect of the 2nd mutation vs. the effect of the first mutation
		vector<vector<double> > temp_stat;
		temp_stat.push_back(y1_t);
		temp_stat.push_back(sd_m);
		
		MultivariateStat ms(temp_stat);
		dir.push_back(ms.regression_slope(1, 0)/mut1_var[t]);
	}
	
	return(dir);	
}

double Direpistasis::individual_fitdir(const DoubleMutantcollection & dmutcol) const
{
	// This is much simpler than for phenotypes, because there is only one dimension. 
	double mut1_var = dmutcol.ref_var_fit();
	// if (mut1_var <= 0.0) cerr << "No variance in fitness among mutants." << endl;
	vector<double> mut1_fit = dmutcol.ref_fit();
	vector<double> mut2_sdfit = dmutcol.var_mutant_fit();
	for (unsigned int i = 0; i < mut2_sdfit.size(); i++) 
	{
		mut2_sdfit[i] = sqrt(mut2_sdfit[i]);
	}
	
	vector<vector<double> > temp_stat;
	temp_stat.push_back(mut1_fit);
	temp_stat.push_back(mut2_sdfit);
	
	MultivariateStat ms(temp_stat);
	return(ms.regression_slope(1, 0)/mut1_var);
}
