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
#include "Population.h"

#include <vector>
#include <string>


/* This class is devoted to the estimation of the canalization (or genetic robustness) of the genetic architecture in a population. 
This is achieved by computing the variance of the effect of random mutations in the population.
The smaller this variance, the more robustness the population displays. 
*/
	
class MiniCanIndiv 
{
	public:
		MiniCanIndiv(const std::vector<Individual>&, const Individual &, bool, bool);
		~MiniCanIndiv() {}
		Phenovec canpheno;
		double canfitness;
};	
	
class Canalization // virtual pure
{
	public:
		Canalization(unsigned int, const Population &, bool logvar = true, bool meancentered = true);
		virtual ~Canalization() = 0;
	
		// get the canalization scores
		Phenovec meanpop_canphen() const;
		Phenovec varpop_canphen() const;		
		double meangene_meanpop_canphen() const;
		std::vector<double> meangene_canphen() const;
		
		double meanpop_canlogfit() const;
		double varpop_canlogfit() const;
		std::vector<double> canlogfit() const;
		
	protected:		
		std::vector<MiniCanIndiv> popcan;		
};
	
class GeneticCanalization : public Canalization
{
	public:
		GeneticCanalization(unsigned int, const Population &, bool logvar = true, bool meancentered = true);
		~GeneticCanalization() { }
};

class DisturbCanalization : public Canalization
{
	public:
		DisturbCanalization(unsigned int, const Population &, bool logvar = true, bool meancentered = true);
		~DisturbCanalization() { }
};

class EnviroCanalization : public Canalization
{
	public:
		EnviroCanalization(unsigned int, const Population &, bool logvar = true, bool meancentered = true);
		~EnviroCanalization() { }
};

#endif // CANALIZATION_H_INCLUDED
