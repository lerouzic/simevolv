// Copyright 2013-2017      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DIREPISTASIS_H_INCLUDED
#define DIREPISTASIS_H_INCLUDED

#include "Population.h"
#include "Individual.h"
#include "Phenotype.h"
#include "Mutantcollection.h"
#include "Fitness.h"

#include <vector>


/* Directional epistasis
	Directional epistasis is estimated by performing two rounds of mutations. 
	1) mutants are generated (the effect of mutation on the phenotype is m1)
	2) previous mutants are mutated again, the variance of the effect of the second mutation is Var(m2|m1) = Var(m1)(1+e m1)^2, where e is the directionality of epistasis. 
	The algorithm performs a linear regression (on the square root of Var(m2|m1) depending on m1), the slope being sqrt(Var(m1)) e m1
*/


class Direpistasis 
{	
	public:
		// The parameters are the total number of mutants, and the population to test
		// Most of the calculation is done in the constructor
		Direpistasis(unsigned int, const Population&);
		~Direpistasis() { }
		
		Phenotype phen_direpistasis() const;
		Phenotype var_phen_direpistasis() const;
		fitness_type fitness_direpistasis() const;
		fitness_type var_fitness_direpistasis() const;
				
	protected:		
		Phenotype individual_direpi(const DoubleMutantcollection &) const;
		fitness_type individual_fitdir(const DoubleMutantcollection &) const;  
		 
		Phenotype dir_epi_mean;
		Phenotype dir_epi_var;
		fitness_type fit_epi_mean;
		fitness_type fit_epi_var;			
};

#endif // DIREPISTASIS_H_INCLUDED
