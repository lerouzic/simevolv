// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Fitness.h"
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;

/* Initialization of the instance (static pointer) */
Fitness * Fitness::instance = NULL;

// constructors and destructor

/* constructor using the parameters from the Pameters files */
/* the Fitness class keeps a (constant) reference of the parameter set.
   (Note: what happens if the parameter set is destroyed? this is not particularly safe) */
Fitness::Fitness(const ParameterSet& param)
    : param(param)
    , type(param.getpar(FITNESS_TYPE)->GetString())
    , strength(param.getpar(FITNESS_STRENGTH)->GetDouble())
    , optimum(param.getpar(FITNESS_OPTIMUM)->GetDouble())
{
}

/* This static function can be used to initialize or reinitialize the Fitness
   single instance with the parameter set.
   Note: I am not sure that the reinitialization from the parameter set can be useful
   in any context. A warning is displayed, but it could be removed if a legitimate use 
   was discovered. */
void Fitness::initialize(const ParameterSet& param)
{
    if (Fitness::instance != NULL)
    {
        delete Fitness::instance;
        Fitness::instance = NULL;
		cerr << "Warning: Fitness is initialized several times." << endl;
    }
    Fitness::instance = new Fitness(param);
}



/* Return the "Population value".
   The meaning of the population value depends on the selection regimes:
   * Linear or exponential selection: the population value is the mean phenotype
   * Truncation selection: the population value is the truncation point
   * Other selection regimes: the population value is meaningless, the function returns 0.0
     but it could be anything. */
     /* Note: this needs to be rethought for multidimensional phenotypes */
double Fitness::GetPopulationValue(const Population& popul)
{
    string type = instance->type;
	if (type == FT_linear || type == FT_expo) 
	{
		return(popul.mean_phenotype());
	} 
	else if (type == FT_truncup) 
	{
		// Broken
		assert("This option is broken");
        //~ sort(pheno.begin(), pheno.end());
        //~ return(pheno[int(instance->strength*pheno.size())]);
     } 
     else if (type == FT_truncdown) 
     {
		 assert("This option is broken");
		//~ sort(pheno.begin(), pheno.end());
        //~ return(pheno[int((1-instance->strength)*pheno.size())]);
     }
     return(0.0);
}


/* This function is called every generation and updates the selection function
   when it is expected to change with time (e.g., fluctuating selection).
   The generation number is used to trigger periodic changes */
void Fitness::update_generation(const long unsigned int generation_number)
{
	// T is the period of fluctuations (for periodic changes)
    long int T = instance->param.getpar(FITNESS_PERIOD)->GetInt();
    
    // Strength of selection. The meaning of s2 depends on the details of the
    // selection regime.
    double s1 = instance->param.getpar(FITNESS_STRENGTH)->GetDouble();
    double s2 = instance->param.getpar(FITNESS_STRENGTH2)->GetDouble();
    
    // Optimum of the selection. The meaning of o2 depends on the details of the
    // selection regime.
    double o1 = instance->param.getpar(FITNESS_OPTIMUM)->GetDouble();
    double o2 = instance->param.getpar(FITNESS_OPTIMUM2)->GetDouble();

	string fluct_type = instance->param.getpar(FITNESS_FLUCT)->GetString();
	if (fluct_type == FF_nofluct) 
	{
		// No fluctuations: the strength and the optimum are directly set from s1 and o1.
		// Note that s2 and o2 have no meaning (reading them from the dataset is useless).
        instance->strength = s1;
        instance->optimum = o1;		
	} 
	else if (fluct_type == FF_smooth) 
	{
		// Smooth fluctuations. The selection optimum follows a sinusoid curve alternating
		// between o1 and o2 with a period T. The selection strength oscillates between s1 and
		// s2 in the same way, whatever it might mean. (be careful to understand the implications
		// of setting different s1 and s2!). 
        instance->strength = s2+(s1-s2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;
        instance->optimum = o2+(o1-o2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;		
	} 
	else if (fluct_type == FF_pflips) 
	{
		// Regular flips between two states (o1 and o2 for the optimum, s1 and s2 for
		// the strength of selection). Every T/2 generations, the selection optimum 
		// shifts from o1 to o2, and T/2 generations afterwards, from o2 to o1.
		// The same for the selection strength, whatever it might mean. 
        if (int(generation_number / double(T/2)) % 2 == 0)
        {
            instance->strength = s1;
            instance->optimum  = o1;
        }
        else
        {
            instance->strength = s2;
            instance->optimum  = o2;
        }		
	} 
	else if (fluct_type == FF_sflips) 
	{
		// The same as for periodic flips, except that o1 / o2 flips (or s1 / s2) occur
		// IN AVERAGE every T generations, but the generation at which it happens is not predictable. 
        // Note that the algorithm is far from perfect: 
        // * 1. it uses a "dirty" trick: every T/2 generations, the strength/optimum is set randomly to 
        //      s1 - o1 or s2 - o2. Thus, half "shifts" do not shift anything, it just remains where it was. 
        //      this allows to write an algorithm that does not need to test its current state. 
        // * 2. when a shift occurs, it occurs for both strength and optimum at the same time. The pairs s1 - o1 and
        //      s2 - o2 are thus never decoupled. I don't know if this is a bug or a feature. 
        if (Random::randnum() < (1.0 / double(T) / 2.0))
        {
            instance->strength = s1;
            instance->optimum  = o1;
        }
        if (Random::randnum() < (1.0 / double(T) / 2.0))
        {
            instance->strength = s2;
            instance->optimum  = o2;
        }		
	} 
	else if (fluct_type == FF_brown) 
	{
		// Brownian motion. Every T generations, the strength and the optimum change, and the change is
		// normally distributed, centered on the former value. 
		// s2 and o2 here stand for the standard deviation of the change. They have a completely different meaning
		// as for other fluctuation regimes. 
        if ((generation_number % T) == 0)
        {
            instance->strength = instance->strength + Random::randgauss()*s2;
            instance->optimum = instance->optimum +  Random::randgauss()*o2;
        }		
	} 
	else if (fluct_type == FF_white) 
	{
		// White noise: the same spirit as for the Brownian motion, but the process has no memory. A new
		// strength and optimum are drawn from a Gaussian distribution every T generations around a constant
		// mean value. 
        if ((generation_number % T) == 0)
        {
            instance->strength = s1 + Random::randgauss()*s2;
            instance->optimum = o1 + Random::randgauss()*o2;
        }		
	} 
	else 
	{
        assert("Fluctuating selection type unknown.");		
	}
}


/* Computes the fitness associated with a given phenotype. A reference to the current population
   has to be provided, because some fitness functions requires a "population value" (see Fitness::GetPopulationValue). 
   This function simply computes this population value, and calls an overloaded version of compute that does
   not need the full population as a reference. */
double Fitness::compute(const Phenotype & phenotype, const Population & population)
{
    double population_value = Fitness::GetPopulationValue(population);
    return(compute(phenotype, population_value));
}


/* Computes the fitness for a particular phenotype, given a "population value", which is necessary for some
   selection regimes. 
   So far, this function is very simple and computes fitness based on only one trait (the first phenotype). This
   will have to change in the close future. */
double Fitness::compute(const Phenotype & phenotype, double population_value)
{
	int focal_phen = 0; // This will have to change for multivariate selection!!!!
	
    assert (instance != NULL);
    double fit = 0.0;

	string type = instance->type;
	if (type == FT_nosel) {
		// No selection: every phenotype has the same fitness (1.0, but it could be any positive number)
		fit = 1.0;
	} else if (type == FT_linear) {
		// Linear selection. w = 1 + s(x-µ), where s is the selection strength, and µ is the population mean. 
		// The population mean is provided as the "population value"
		// Note that linear selection could lead to negative values of fitness, which need to be set to 0.
        fit = 1.0 + instance->strength*(phenotype[focal_phen] - population_value);
        if (fit < 0)
        {
            fit = 0;
        }		
	} else if (type == FT_expo) {
		// Exponential selection. w = exp[s(x-µ)]
		fit = exp(instance->strength*(phenotype[focal_phen] - population_value));
	} else if (type == FT_gauss) {
		// Gaussian selection (stabilizing around an optimum o).
		// w = exp[-s(x-o)²]
		fit = exp(-instance->strength*
			(phenotype[focal_phen] - instance->optimum)*
            (phenotype[focal_phen] - instance->optimum));		
	} else if (type == FT_quad) {
		// Quadratic selection: w = 1 - s(x-o)².
		// This can lead to negative fitness values that must be set to 0. 
        fit = 1.0 - instance->strength*
              (phenotype[focal_phen] - instance->optimum)*
              (phenotype[focal_phen] - instance->optimum);
        if (fit < 0)
        {
            fit = 0;
        }		
	} else if (type == FT_truncup) {
		// Up truncation selection. w = 1 if x > t, w = 0 otherwise.
		// The threshold t is provided as the "population value".
        fit = 0.0;
        if (phenotype[focal_phen] > population_value)
            fit = 1.0;		
	} else if (type == FT_truncdown) {
		// Down truncation selection. w = 1 if w < t, w = 0 otherwise. 
        fit = 0.0;
        if (phenotype[focal_phen] < population_value)
            fit = 1.0;		
	} else if (type == FT_concave) {
		// Concave fitness function. w = 1 + log(1 + 2s(x-µ))/2
		// The slope of this function is s when x = µ
		// Note that this function can return negative fitnesses. 
        fit = 1.0 + 0.5*log(1.0+2.0*instance->strength*(phenotype[focal_phen] - population_value));
        if (fit < 0.0)
            fit = 0.0;		
	} else if (type == FT_convex) {
		// Convex fitness function. w = exp(-sqrt(s)|x-µ|)
		// As the previous one, the slope is s when x = µ. 
        fit = exp(-sqrt(instance->strength*
            (phenotype[focal_phen] - instance->optimum)*
            (phenotype[focal_phen] - instance->optimum)));		
	} else {
		assert("Fitness type unknown.");
	}

    return(fit);
}

