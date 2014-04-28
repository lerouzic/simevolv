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


#include "Phenotype.h" 

#include <vector>



class Population;


class Heritability 
{
	public:
		// constructors
		Heritability(unsigned int, const Population &);
		
		// functions
		Phenotype h2() const;
		Phenotype H2() const;
		double fit_h2() const;
		
	protected:
		struct ParentOffspring 
		{
			Phenotype mother_phen;
			Phenotype mother_gen;
			double mother_fit;
			Phenotype father_phen;
			Phenotype father_gen;
			double father_fit;
			Phenotype offspring_phen;
			Phenotype offspring_gen;
			double offspring_fit;
		};
		std::vector<ParentOffspring> parentoffspring;	
};

#endif // HERITABILITY_H_INCLUDED
