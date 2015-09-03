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
		double fitness;
};


class Mutantcollection
{
	friend class DoubleMutantcollection;
	
	public: 
		Mutantcollection(unsigned int, const Individual &, const Population &);
		~Mutantcollection();
		
		Mutantcollection & operator=(const Mutantcollection &);
		
		Phenovec mean_phen() const;
		Phenovec var_phen() const;
		double mean_fit() const;
		double var_fit() const;
		
	protected:
		void compute_phenostat() const;
		void compute_fitstat() const;
	
		MiniIndividual reference;
		std::vector<MiniIndividual> collection;
		
		mutable PhenotypeStat* phenostat;
		mutable UnivariateStat* fitstat;
};
	
	
class DoubleMutantcollection
{
	public:
		DoubleMutantcollection(unsigned int, unsigned int, const Individual &, const Population &);
		~DoubleMutantcollection();
		
		DoubleMutantcollection & operator=(const DoubleMutantcollection &);
		
		Phenovec ref_mean_phen() const;
		Phenovec ref_var_phen() const;
		double ref_mean_fit() const;
		double ref_var_fit() const;
		
		std::vector<Phenovec> ref_phen() const;
		std::vector<Phenovec> var_mutant_phen() const;		
		std::vector<double> ref_fit() const;
		std::vector<double> var_mutant_fit() const;
		
	protected:
		void compute_refstat() const;	
		void compute_reffitstat() const;
					
		std::vector<Mutantcollection> dcollection;
		
		mutable PhenotypeStat* refstat;
		mutable UnivariateStat* reffitstat;
};	
	
#endif // MUTANTCOLLECTION_H_INCLUDED
