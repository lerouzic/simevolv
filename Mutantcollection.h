// Copyright 2013-2014      Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef MUTANTCOLLECTION_H_INCLUDED
#define MUTANTCOLLECTION_H_INCLUDED

#include "Phenotype.h"
#include "Population.h"
#include "Individual.h"
#include "Statistics.h"

#include <vector>


struct MiniIndividual
{
	public:
		MiniIndividual(const MiniIndividual & mini) : phen(mini.phen), fitness(mini.fitness) { }
		MiniIndividual(const Individual & ind) : phen(ind.get_phenotype()), fitness(ind.get_fitness()) { }
				
		Phenotype phen;
		basic_fitness fitness;
};


class Mutantcollection
{
	friend class DoubleMutantcollection;
	
	public: 
		Mutantcollection(unsigned int, const Individual &, const Population &);
		~Mutantcollection();
		
		Mutantcollection & operator=(const Mutantcollection &);
		
		Phenotype mean_phen() const;
		Phenotype var_phen() const;
		basic_fitness mean_fit() const;
		basic_fitness var_fit() const;
		
	protected:
		MiniIndividual reference;
		std::vector<MiniIndividual> collection;
};
	
	
class DoubleMutantcollection
{
	public:
		DoubleMutantcollection(unsigned int, unsigned int, const Individual &, const Population &);
		~DoubleMutantcollection();
		
		DoubleMutantcollection & operator=(const DoubleMutantcollection &);
		
		Phenotype ref_mean_phen() const;
		Phenotype ref_var_phen() const;
		basic_fitness ref_mean_fit() const;
		basic_fitness ref_var_fit() const;
		
		std::vector<Phenotype> ref_phen() const;
		std::vector<Phenotype> var_mutant_phen() const;		
		std::vector<basic_fitness> ref_fit() const;
		std::vector<basic_fitness> var_mutant_fit() const;
		
	protected:
		std::vector<Mutantcollection> dcollection;
};	
	
#endif // MUTANTCOLLECTION_H_INCLUDED
