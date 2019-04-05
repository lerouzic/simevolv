// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef HERITABILITY_H_INCLUDED
#define HERITABILITY_H_INCLUDED

#include "types.h"
#include "Phenotype.h" 
#include "Population.h"
#include "Fitness.h"

#include <vector>
 
class Heritability 
{
	public:
		// constructors
		Heritability(unsigned int, const Population &); // parameters are the number of parent/offspring pairs to be tested, and the population
		~Heritability() { }
		
		// functions
		Phenotype h2() const;	// Narrow-sense heritability (VA/VP)
		fitness_type fit_h2() const;	// Narrow-sense heritability in fitness
		Phenotype covOffPar(size_t trait) const; 	// Covariance between all offspring traits and a mid-parent trait
		
	protected:
		// internal structure for parent - offspring pairs (no need to make this visible from outside)
		struct ParentOffspring 
		{
			Phenotype mother_phen;
			Phenotype mother_gen;
			fitness_type mother_fit;
			Phenotype father_phen;
			Phenotype father_gen;
			fitness_type father_fit;
			Phenotype offspring_phen;
			Phenotype offspring_gen;
			fitness_type offspring_fit;
		};
		std::vector<ParentOffspring> parentoffspring;	
};

#endif // HERITABILITY_H_INCLUDED
