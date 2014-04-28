// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef CANALIZATION_H_INCLUDED
#define CANALIZATION_H_INCLUDED

#include "Phenotype.h"
#include "Statistics.h"
#include "Parameters.h"

#include <vector>



class Individual; // crossed references: Population.h -> Canalization.h -> Individual.h -> Population.h
class Population;


class Canalization
{
	public:
		Canalization(unsigned int, const Population &);
						
		// get the results
		Phenotype phen_canalization();
		double fitness_canalization();
				
	protected:
		// fill the object
		void reference_indiv(Individual);
		void mutant_indiv (Individual);
		
		void process();	
		void process_phen();
		void process_fit();
		
		std::vector<Individual> reference;
		std::vector<std::vector<Individual> > mutants;
		
		std::vector<Phenotype> mean_per_indiv;
		Phenotype mean_of_mean;
		Phenotype var_of_mean;
		
		std::vector<Phenotype> var_per_indiv;
		Phenotype mean_of_var;
		Phenotype var_of_var;
		
		std::vector<double> indiv_fitness_mean;
		std::vector<double> indiv_fitness_var;
		
		bool phen_ready;
		bool fit_ready;
};


#endif // CANALIZATION_H_INCLUDED
