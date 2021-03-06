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
#include <gsl/gsl_permutation.h> // for matrix inversion
#include <gsl/gsl_linalg.h>

using namespace std;

/* Initialization of the instance (static pointer) */
Fitness * Fitness::instance = NULL;

/* constructor using the parameters from the Pameters files */
Fitness::Fitness(const ParameterSet& param)
{
	string type = param.getpar(FITNESS_TYPE)->GetString();
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
	else if (type == FT_multigauss)
	{
		fitphen = new Fitness_Phenotype_MultivarGaussian(param);
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

void Fitness::Terminate() 
{
    if (instance) 
        delete instance;
    instance = NULL;
}

/*  The destructor does something here, as it cleans all pointers towards fitness classes 
 * (Question: is this destructor called at all?) */
Fitness::~Fitness()
{
	delete fitphen;
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


Fitness_Phenotype_MultivarGaussian::Fitness_Phenotype_MultivarGaussian(const ParameterSet & param)
	: Fitness_Phenotype_Stabilizing(param)
	, cor(param.getpar(FITNESS_CORRELATION)->GetVectorDouble())
	, invsigma(nullptr)
{
}

Fitness_Phenotype_MultivarGaussian::~Fitness_Phenotype_MultivarGaussian() 
{
	gsl_matrix_free(invsigma);
	invsigma = nullptr;
}


fitness_type Fitness_Phenotype_MultivarGaussian::get_fitness(const Phenotype& phenotype, const Population& population) {
	assert(phenotype.is_defined());
	
	size_t n_traits = phenotype.dimensionality();
	
	// Check if all vectors have been properly extended (first call only).
	if (strength.size() < n_traits) strength = expand_vec(strength, n_traits);
	if (optimum.size() < n_traits) optimum = expand_vec(optimum, n_traits);
	if (cor.size() < ((n_traits*n_traits) - n_traits)/2) cor = expand_vec(cor, (n_traits*n_traits - n_traits)/2);
	
	// Using gsl_matrix, as gsl is already a dependency of the program. This algorithm thus uses C-style code.
	
	if (invsigma == nullptr) { // This should be run only once (or each time something changes in the selection)
		compute_invsigma(n_traits);
	}
	
	fitness_type fit = 0.0;
	
	for (size_t i = 0; i < n_traits; i++) {
		for (size_t j = 0; j < n_traits; j++) {
			pheno_type diff1 = phenotype[i] - optimum[i];
			pheno_type diff2 = phenotype[j] - optimum[j];
			fit += diff1 * diff2 * gsl_matrix_get(invsigma, i, j);
		}
	}
	return(exp(-fit/2.0));
}

void Fitness_Phenotype_MultivarGaussian::compute_invsigma(const size_t n_traits) {
	if (strength.size() < n_traits) strength = expand_vec(strength, n_traits);
	if (cor.size() < ((n_traits^2) - n_traits)/2) cor = expand_vec(cor, (n_traits*n_traits - n_traits)/2);
	
	// when inverting the selection matrix, zero strengths -> inf variances and the matrix is not positive definite. 
	for (auto &i : strength) if (i == 0.0) i = 1e-20; // numeric_limits<fitness_type>::min();
	
	// The first step is to reconstruct a variance-covariance matrix from selection strengths and correlations
	gsl_matrix *mat = gsl_matrix_alloc(n_traits, n_traits);
	
	for (size_t i = 0; i < n_traits; i++) {
		gsl_matrix_set(mat, i, i, 1./2./abs(strength[i]));
		if (i < n_traits - 1) {
			for (size_t j = i+1; j < n_traits; j++) {
				size_t cor_index = i+j*(j-1)/2; // index of element i, j in a upper triangular matrix
				fitness_type cov = cor[cor_index]*sqrt(1./2./abs(strength[i]))*sqrt(1./2./abs(strength[j]));
				gsl_matrix_set(mat, i, j, cov);
				gsl_matrix_set(mat, j, i, cov);
			}
		}
	}
	
	// We then need to invert this matrix
	gsl_permutation *p = gsl_permutation_alloc(n_traits);
	int s;
	gsl_linalg_LU_decomp(mat, p, &s); // Hoping that the way we build the matrix ensures positive definite propreties. Otherwise, crash. 
	invsigma = gsl_matrix_alloc(n_traits, n_traits);
	gsl_linalg_LU_invert(mat, p, invsigma);
	gsl_permutation_free(p);
	gsl_matrix_free(mat);
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
