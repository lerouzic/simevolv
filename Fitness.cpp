// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2017 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Fitness.h"

#include "Population.h"
#include "Parconst.h"
#include "Random.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

/* Initialization of the instance (static pointer) */
Fitness * Fitness::instance = NULL;

/* constructor using the parameters from the Pameters files */
Fitness::Fitness(const ParameterSet& param)
{
	string type = param.getpar(FITNESS_TYPE)->GetString();
	string fluct_type = param.getpar(FITNESS_FLUCT)->GetString();
	string stab_type = param.getpar(FITNESS_STAB)->GetString();
		
	if (type == FT_nosel) 
	{
		fitphen = new Fitness_Phenotype_Noselection(param);
	} 
	else if (type == FT_linear) 
	{
		fitphen = new Fitness_Phenotype_Linear(param);
	} 
	else if (type == FT_expo) 
	{
		fitphen = new Fitness_Phenotype_Expo(param);
	} 
	else if (type == FT_concave) 
	{
		fitphen = new Fitness_Phenotype_Concave(param);
	} 
	else if (type == FT_gauss) 
	{
		fitphen = new Fitness_Phenotype_Gaussian(param);
	} 
	else if (type == FT_quad) 
	{
		fitphen = new Fitness_Phenotype_Quadratic(param);
	} 
	else if (type == FT_convex) 
	{
		fitphen = new Fitness_Phenotype_Biconvex(param);
	} 
	else 
	{
		assert("Fitness type unknown or not implemented");
	}
	
	if (fluct_type == FF_nofluct) 
	{
		fitfluct = new Fitness_Fluct_Nofluct(param);
	} 
	else if (fluct_type == FF_smooth) 
	{
		fitfluct = new Fitness_Fluct_Smooth(param);
	} 
	else if (fluct_type == FF_pflips) 
	{
		fitfluct = new Fitness_Fluct_Flips(param);
	} 
	else if (fluct_type == FF_sflips) 
	{
		fitfluct = new Fitness_Fluct_Sflips(param);
	} 
	else if (fluct_type == FF_brown) 
	{
		fitfluct = new Fitness_Fluct_Brownian(param);
	} 
	else if (fluct_type == FF_white) 
	{
		fitfluct = new Fitness_Fluct_Whitenoise(param);
	} 
	else 
	{
		assert("Fluctuation type unknown or not implemented");
	}
	
	if (stab_type == FS_nostab) 
	{
		fitstab = new Fitness_Stability_Noselection(param);
	} 
	else if (stab_type == FS_expo) 
	{
		fitstab = new Fitness_Stability_Exponential(param);
	} 
	else 
	{
		assert("Stability selection type unknown or not implemented");
	}
}

/*  The destructor does something here, as it cleans all pointers towards fitness classes 
 * (Question: is this destructor called at all?) */
Fitness::~Fitness()
{
	delete fitphen;
	delete fitfluct;
	delete fitstab;
}

/* This static function can be used to initialize or reinitialize the Fitness
   single instance with the parameter set. */
void Fitness::initialize(const ParameterSet& param)
{
    if (Fitness::instance != NULL)
    {
        delete Fitness::instance;
        Fitness::instance = NULL;
		// cerr << "Warning: Fitness is initialized several times." << endl;
    }
    Fitness::instance = new Fitness(param);
}

/* This function should be called every generation, before computing individual fitnesses.
   It ensures that the fitness class is initialized with the current population, which
   is necessary for some fitness types. Note that there is an (imperfect!) control if
   update is not called prior to computing fitnesses, but one should never rely on it! */
void Fitness::update(const Population & population) 
{
	Fitness::instance->fitphen->update(population);
}

/* This function is called every generation and updates the selection function
   when it is expected to change with time (e.g., fluctuating selection).
   The generation number is used to trigger periodic changes */
void Fitness::fluctuate(unsigned int generation_number)
{
	Fitness::instance->fitphen->fluctuate(Fitness::instance->fitfluct, generation_number);
}

/* Computes the fitness associated with a given phenotype. A reference to the current population
   has to be provided, because some fitness functions requires a "population value".
   Fitness is the product of two components: a fitness score from the phenotypic values, 
   and a fitness score on stability (which has a meaning only for some genetic architectures). */
fitness_type Fitness::compute(const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	return(Fitness::instance->fitphen->get_fitness(phenotype, population) 
		* Fitness::instance->fitstab->get_fitness(phenotype));
}

FitnessOptimum Fitness::current_optimum()
{
	return(Fitness::instance->fitphen->get_optimum());
}


/******************************** Fitness Fluct ******************************/

Fitness_Fluct::~Fitness_Fluct()
{ // virtual pure
}


//////// Fitness_Fluct_Nofluct

Fitness_Fluct_Nofluct::Fitness_Fluct_Nofluct(const ParameterSet & param)
	: Fitness_Fluct(param)
{
}

FitnessStrength Fitness_Fluct_Nofluct::get_new_strength(const FitnessStrength & old_strength, unsigned int generation) 
{ 
	return(old_strength); 
}

FitnessOptimum Fitness_Fluct_Nofluct::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation) 
{ 
	return(old_optimum); 
} 


//////// Fitness_Fluct_States

Fitness_Fluct_States::~Fitness_Fluct_States()
{ // The program does not compile without the empty virtual pure destructor
}

Fitness_Fluct_States::Fitness_Fluct_States(const ParameterSet & param)
	: Fitness_Fluct(param)
	, strength_state1(param.getpar(FITNESS_STRENGTH)->GetVectorDouble())
	, strength_state2(param.getpar(FITNESS_STRENGTH2)->GetVectorDouble())
	, optimum_state1(param.getpar(FITNESS_OPTIMUM)->GetVectorDouble())
	, optimum_state2(param.getpar(FITNESS_OPTIMUM2)->GetVectorDouble())
	, period(param.getpar(FITNESS_PERIOD)->GetInt())
{
}

Fitness_Fluct_Flips::Fitness_Fluct_Flips(const ParameterSet & param)
	: Fitness_Fluct_States(param)
{
}

FitnessStrength Fitness_Fluct_Flips::get_new_strength(const FitnessStrength & old_strength, unsigned int generation)
{ // old_strength is used only to compare vector lengths
	if (old_strength.size() > strength_state1.size()) 
	{
		strength_state1 = expand_vec(strength_state1, old_strength.size());
		strength_state2 = expand_vec(strength_state2, old_strength.size());
	}
	if (int(generation / float(period/2)) % 2 == 0)
	{
        return(strength_state1);
    } 
    else 
    {
        return(strength_state2);
    }
}

FitnessOptimum Fitness_Fluct_Flips::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation)
{ // old_optimum is used only to compare vector lengths
	if (old_optimum.size() > optimum_state1.size()) 
	{
		optimum_state1 = expand_vec(optimum_state1, old_optimum.size());
		optimum_state2 = expand_vec(optimum_state2, old_optimum.size());
	}
	if (int(generation / float(period/2)) % 2 == 0)
	{
        return(optimum_state1);
    } 
    else 
    {
        return(optimum_state2);
    }
}

Fitness_Fluct_Sflips::Fitness_Fluct_Sflips(const ParameterSet & param)
	: Fitness_Fluct_States(param)
{	
}

FitnessStrength Fitness_Fluct_Sflips::get_new_strength(const FitnessStrength & old_strength, unsigned int generation)
{ // old_strength is used only to compare vector lengths
	if (old_strength.size() > strength_state1.size()) 
	{
		strength_state1 = expand_vec(strength_state1, old_strength.size());
		strength_state2 = expand_vec(strength_state2, old_strength.size());
	}

	if (Random::randnum() < (1.0 / float(period) / 2.0))
	{
        return(strength_state1);
    } 
    else 
    {
        return(strength_state2);
    }
}

FitnessOptimum Fitness_Fluct_Sflips::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation)
{ // old_optimum is used only to compare vector lengths
	if (old_optimum.size() > optimum_state1.size()) 
	{
		optimum_state1 = expand_vec(optimum_state1, old_optimum.size());
		optimum_state2 = expand_vec(optimum_state2, old_optimum.size());
	}
	
	if (Random::randnum() < (1.0 / float(period) / 2.0))
	{
        return(optimum_state1);
    } 
    else 
    {
        return(optimum_state2);
    }
}

Fitness_Fluct_Smooth::Fitness_Fluct_Smooth(const ParameterSet & param)
	: Fitness_Fluct_States(param)
{
}

FitnessStrength Fitness_Fluct_Smooth::get_new_strength(const FitnessStrength & old_strength, unsigned int generation)
{
	if (old_strength.size() > strength_state1.size()) 
	{
		strength_state1 = expand_vec(strength_state1, old_strength.size());
		strength_state2 = expand_vec(strength_state2, old_strength.size());
	}
	FitnessStrength new_strength;
	for (unsigned int tt = 0; tt < old_strength.size(); tt++) 
	{
		new_strength.push_back(strength_state2[tt]+(strength_state1[tt]-strength_state2[tt])*(1.0+cos(2.0*generation*M_PI/float(period)))/2.0);
	}
	return(new_strength);
}

FitnessOptimum Fitness_Fluct_Smooth::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation)
{
	if (old_optimum.size() > optimum_state1.size()) 
	{
		optimum_state1 = expand_vec(optimum_state1, old_optimum.size());
		optimum_state2 = expand_vec(optimum_state2, old_optimum.size());
	}
	FitnessOptimum new_optimum;
	for (unsigned int tt = 0; tt < old_optimum.size(); tt++) 
	{
		new_optimum.push_back(optimum_state2[tt]+(optimum_state1[tt]-optimum_state2[tt])*(1.0+cos(2.0*generation*M_PI/float(period)))/2.0);
	}
	return(new_optimum);
}


//////// Fitness_Fluct_Noise

Fitness_Fluct_Noise::Fitness_Fluct_Noise(const ParameterSet & param)
	: Fitness_Fluct(param)
	, strength_sd(param.getpar(FITNESS_STRENGTH2)->GetVectorDouble())
	, optimum_sd(param.getpar(FITNESS_OPTIMUM2)->GetVectorDouble())
	, period(param.getpar(FITNESS_PERIOD)->GetInt())
{	
}

Fitness_Fluct_Noise::~Fitness_Fluct_Noise()
{
}

Fitness_Fluct_Whitenoise::Fitness_Fluct_Whitenoise(const ParameterSet & param)
	: Fitness_Fluct_Noise(param)
	, strength_ref(param.getpar(FITNESS_STRENGTH)->GetVectorDouble())
	, optimum_ref(param.getpar(FITNESS_OPTIMUM)->GetVectorDouble())
{
}

FitnessStrength Fitness_Fluct_Whitenoise::get_new_strength(const FitnessStrength & old_strength, unsigned int generation)
{
	if (old_strength.size() > strength_ref.size()) 
	{
		strength_ref = expand_vec(strength_ref, old_strength.size());
		strength_sd = expand_vec(strength_sd, old_strength.size());
	}
	if (generation % period == 0) 
	{
		FitnessStrength new_strength;
		for (unsigned int tt = 0; tt < old_strength.size(); tt++) 
		{
			new_strength.push_back(strength_ref[tt] + strength_sd[tt]*Random::randgauss());
		}
		return(new_strength);
	} 
	else 
	{
		return(old_strength);
	}
}

FitnessOptimum Fitness_Fluct_Whitenoise::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation)
{
	if (old_optimum.size() > optimum_ref.size()) 
	{
		optimum_ref = expand_vec(optimum_ref, old_optimum.size());
		optimum_sd = expand_vec(optimum_sd, old_optimum.size());
	}
	if (generation % period == 0) 
	{
		FitnessOptimum new_optimum;
		for (unsigned int tt = 0; tt < old_optimum.size(); tt++) 
		{
			new_optimum.push_back(optimum_ref[tt] + optimum_sd[tt]*Random::randgauss());
		}
		return(new_optimum);
	} 
	else 
	{
		return(old_optimum);
	}
}

Fitness_Fluct_Brownian::Fitness_Fluct_Brownian(const ParameterSet & param)
	: Fitness_Fluct_Noise(param)
{	
}

FitnessStrength Fitness_Fluct_Brownian::get_new_strength(const FitnessStrength & old_strength, unsigned int generation)
{
	if (old_strength.size() > strength_sd.size()) 
	{
		strength_sd = expand_vec(strength_sd, old_strength.size());
	}
	if (generation % period == 0) 
	{
		FitnessStrength new_strength;
		for (unsigned int tt = 0; tt < old_strength.size(); tt++) 
		{
			new_strength.push_back(old_strength[tt] + strength_sd[tt]*Random::randgauss());
		}
		return(new_strength);
	} 
	else 
	{
		return(old_strength);
	}
}

FitnessOptimum Fitness_Fluct_Brownian::get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation)
{
	if (old_optimum.size() > optimum_sd.size()) 
	{
		optimum_sd = expand_vec(optimum_sd, old_optimum.size());
	}
	if (generation % period == 0) 
	{
		FitnessOptimum new_optimum;
		for (unsigned int tt = 0; tt < old_optimum.size(); tt++) 
		{
			new_optimum.push_back(old_optimum[tt] + optimum_sd[tt]*Random::randgauss());
		}
		return(new_optimum);
	} 
	else 
	{
		return(old_optimum);
	}
}



/******************************** Fitness_Phenotype ***************************/

Fitness_Phenotype::~Fitness_Phenotype() 
{ // Empty virtual pure destructor
	
}

fitness_type Fitness_Phenotype::get_fitness(const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	unsigned int dim = phenotype.dimensionality();
	fitness_type fitness = 1.0;
	for (unsigned int dd = 0; dd < dim; dd++) 
	{
		fitness *= get_fitness_trait(dd, phenotype, population);
	}
	return(fitness);
}

FitnessOptimum Fitness_Phenotype::get_optimum() const 
{ 
	FitnessOptimum nothing; 
	return(nothing); 
}


//////// Fitness_Phenotype_Noselection

Fitness_Phenotype_Noselection::Fitness_Phenotype_Noselection(const ParameterSet & param) 
	: Fitness_Phenotype(param)
{	
}


//////// Fitness_Phenotype_Directional

Fitness_Phenotype_Directional::Fitness_Phenotype_Directional(const ParameterSet & param)
	: Fitness_Phenotype(param)
	, strength(param.getpar(FITNESS_STRENGTH)->GetVectorDouble())
	, popmem(NULL)
{
}

Fitness_Phenotype_Directional::~Fitness_Phenotype_Directional()
{ // virtual pure
}

void Fitness_Phenotype_Directional::update(const Population & population)
{
	popmean = population.mean_phenotype();
	popmem = &population;
}

void Fitness_Phenotype_Directional::fluctuate(Fitness_Fluct * fitfluct, unsigned int generation)
{
	strength = fitfluct->get_new_strength(strength, generation);
}

Fitness_Phenotype_Linear::Fitness_Phenotype_Linear(const ParameterSet & param)
	: Fitness_Phenotype_Directional(param)
{	
}

fitness_type Fitness_Phenotype_Linear::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
		
	if (popmem != &population) 
	{
		cerr << "Warning, Fitness was not updated after a population change." << endl;
		update(population);
	}
	fitness_type fit = strength[trait]*(phenotype[trait] - popmean[trait]);
	if (fit < 0.0)
	{ 
		fit = 0.0;
	}
	return(fit); 
}

Fitness_Phenotype_Expo::Fitness_Phenotype_Expo(const ParameterSet & param)
	: Fitness_Phenotype_Directional(param)
{
}

fitness_type Fitness_Phenotype_Expo::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
		
	if (popmem != &population) 
	{
		cerr << "Warning, Fitness was not updated after a population change." << endl;
		update(population);
	}
	pheno_type departure = phenotype[trait] - popmean[trait];
	return(exp(strength[trait]*departure)); 
}

Fitness_Phenotype_Concave::Fitness_Phenotype_Concave(const ParameterSet & param)
	: Fitness_Phenotype_Directional(param)
{
}

fitness_type Fitness_Phenotype_Concave::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
		
	if (popmem != &population) 
	{
		cerr << "Warning, Fitness was not updated after a population change." << endl;
		update(population);
	}
	fitness_type fit = 1.0 + 0.5*log(1.0+2.0*strength[trait]*(phenotype[trait] - popmean[trait]));
	if (fit < 0.0) 
	{
		fit = 0.0;
	}
	return(fit); 
}


//////// Fitness_Phentoype_Stabilizing

Fitness_Phenotype_Stabilizing::Fitness_Phenotype_Stabilizing(const ParameterSet & param) 
	: Fitness_Phenotype(param)
	, strength(param.getpar(FITNESS_STRENGTH)->GetVectorDouble())
	, optimum(param.getpar(FITNESS_OPTIMUM)->GetVectorDouble())
{
}

void Fitness_Phenotype_Stabilizing::fluctuate(Fitness_Fluct * fitfluct, unsigned int generation)
{
	strength = fitfluct->get_new_strength(strength, generation);
	optimum = fitfluct->get_new_optimum(optimum, generation);
}

Fitness_Phenotype_Gaussian::Fitness_Phenotype_Gaussian(const ParameterSet & param)
	: Fitness_Phenotype_Stabilizing(param)
{
}

fitness_type Fitness_Phenotype_Gaussian::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	assert(trait < phenotype.dimensionality());
	
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
	
	if (optimum.size() < phenotype.dimensionality())
	{
		optimum = expand_vec(optimum, phenotype.dimensionality());
        
	}
	
	pheno_type departure = phenotype[trait] - optimum[trait];
	fitness_type fit = exp(- strength[trait]*departure*departure);
	return(fit);
}

Fitness_Phenotype_Quadratic::Fitness_Phenotype_Quadratic(const ParameterSet & param) 
	: Fitness_Phenotype_Stabilizing(param)
{
}

fitness_type Fitness_Phenotype_Quadratic::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
    assert(phenotype.is_defined());
	assert(trait < phenotype.dimensionality());
	
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
	
	if (optimum.size() < phenotype.dimensionality())
	{
		optimum = expand_vec(optimum, phenotype.dimensionality());	
	}

	pheno_type departure = phenotype[trait] - optimum[trait];
	fitness_type fit = 1.0 - strength[trait]*departure*departure;
	if (fit < 0.0)
	{
		fit = 0.0;
	}
	
	return(fit);
} 

Fitness_Phenotype_Biconvex::Fitness_Phenotype_Biconvex(const ParameterSet & param)
	: Fitness_Phenotype_Stabilizing(param)
{
}

fitness_type Fitness_Phenotype_Biconvex::get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population & population)
{
	assert(trait < phenotype.dimensionality());
	
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
	
	if (optimum.size() < phenotype.dimensionality())
	{
		optimum = expand_vec(optimum, phenotype.dimensionality());	
	}

	pheno_type departure = phenotype[trait] - optimum[trait];
	fitness_type fit = exp(-sqrt(strength[trait]*strength[trait]*departure*departure));
	return(fit);	
}



/********************** Fitness_Stability *******************/

fitness_type Fitness_Stability::get_fitness(const Phenotype & phenotype)
{
	unsigned int dim = phenotype.dimensionality();
	fitness_type answer = 1.0;
	for (unsigned int dd = 0; dd < dim; dd++) 
	{
		answer *= get_fitness_trait(dd, phenotype);
	} 
	return(answer);
}

Fitness_Stability_Noselection::Fitness_Stability_Noselection(const ParameterSet & param)
	: Fitness_Stability(param)
{
}

fitness_type Fitness_Stability_Noselection::get_fitness_trait(unsigned int trait, const Phenotype & phenotype) 
{ 
	return(1.0); 
}

Fitness_Stability_Exponential::Fitness_Stability_Exponential(const ParameterSet & param) 
	: Fitness_Stability(param)
	, strength(param.getpar(FITNESS_STABSTR)->GetVectorDouble())
{	
}

fitness_type Fitness_Stability_Exponential::get_fitness_trait(unsigned int trait, const Phenotype & phenotype)
{
	assert(trait < phenotype.dimensionality());
	
	if (strength.size() < phenotype.dimensionality())
	{
		strength = expand_vec(strength, phenotype.dimensionality());
	}
		
	return(exp(-strength[trait]*phenotype.get_pheno2(trait)));
}
