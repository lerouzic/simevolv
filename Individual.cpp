// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Individual.h"

#include "Population.h"
#include "Fitness.h"
#include "Architecture.h"
#include "Parconst.h"

#include <vector>
#include <sstream>

using namespace std;



// constructors and destructor

/* constructor using two haplotypes */
Individual::Individual(const Haplotype& gam_father, const Haplotype& gam_mother)
    : genotype(gam_father, gam_mother)
{
    initialize();
}

Individual::Individual(const Haplotype& gam_father, const Haplotype& gam_mother, 
	const EpigeneticInfo& epimother)
	: genotype(gam_father, gam_mother), epiinfo(epimother)
{
	initialize();
}

/* constructor using the parameters from ParameterSet */
Individual::Individual(const ParameterSet& param)
	: genotype(param)
	, epigenet(param.getpar(GENET_EPIGENET)-> GetDouble())
{
    initialize();
}

Individual::Individual(const Individual& copy)
    : genotype(copy.genotype)
    , phenotype(copy.phenotype)
    , fitness(copy.fitness)
    , epigenet(copy.epigenet)
    , epiinfo(copy.epiinfo)
{
}

Individual::~Individual()
{
}

// operator overload

Individual & Individual::operator = (const Individual& copy)
{
    if (this == &copy)
        return (*this);

    genotype=copy.genotype;
    phenotype=copy.phenotype;
    fitness=copy.fitness;
    epigenet=copy.epigenet;
    epiinfo=copy.epiinfo;

    return(*this);
}


// instance and initialization

/* initialize the individual, with genotypic, phenotypic and environmental values */
void Individual::initialize()
{
    Architecture * archi = Architecture::Get();
    phenotype = archi -> phenotypic_value(genotype, true, epiinfo);
    fitness = 0;
    
    // So far, the epigenetic factor does not evolve).
    // If epiinfo is not initialized (no mother = first generation)
    // its value has already been set from the parameter set. 
    // otherwise, we just copy the mother's factor. 
    if (epiinfo.is_defined())
		epigenet = epiinfo.get_epigenet();
}

// functions
void Individual::update_fitness(const Population & pop)
{
    fitness = Fitness::compute(phenotype, pop);
}


// getters
double Individual::get_fitness() const
{
    return(fitness);
}

double Individual::get_epigenet() const
{
	return(epigenet);
}

Phenotype Individual::get_phenotype() const
{
    return(phenotype);
}

string Individual::write_debug(unsigned int gam) const 
{
	assert ((gam==1) || (gam==2));
	// Dirty but convenient: write any debug information into a string.
	ostringstream o;
	
	if (gam==1) 
		o << genotype.gam_father.write_debug();
	else
		o << genotype.gam_mother.write_debug();
	
	return(o.str());
}

EpigeneticInfo Individual::make_epiinfo() const
// Object EpigeneticInfo created by the mother
// to be stored in the offspring
{
	EpigeneticInfo epi(epigenet, phenotype);
	return(epi);
}

/* create a new individual from the paternal and maternal gametes */
Individual Individual::mate(const Individual& father, const Individual& mother)
{
    Individual offspring(father.produce_gamete(), mother.produce_gamete(), mother.make_epiinfo());
    return(offspring);
}

/* produce the gametes of an individual : recombination and mutation */
Haplotype Individual::produce_gamete() const
{
    Haplotype gamete(genotype.recombine());
    gamete.draw_mutation();
    return(gamete);
}

/* determines if there will be a mutation in the individual */
/* This should be called in special cases (like when generating mutant clones)
 * but not for regular simulations during which draw_mutation() is called on the
 * gametes */
/* Note: mutations on individuals change the genotype and the genotypic value.
 * The phenotype is updated assuming no environmental noise (?!?)
 * The fitness is set to 0 */
void Individual::draw_mutation()
{
    genotype.draw_mutation();
    phenotype = Architecture::Get() -> phenotypic_value(genotype, true, epiinfo);
    fitness = 0.0;
}

/* force to make a mutation in an individual */
void Individual::make_mutation(bool test /* = false */)
{
    genotype.make_mutation(test);
    phenotype = Architecture::Get() -> phenotypic_value(genotype, true, epiinfo);
    fitness = 0.0;
} 

/* test the canalization : produce the clones used for the calculation */
Individual Individual::test_canalization(unsigned int nb_mut, const Population & pop) const 
{
	Individual clone(*this);
	for (unsigned int mut = 0; mut < nb_mut; mut++) 
	{
		clone.make_mutation(true);
	}
	clone.update_fitness(pop);	
	return(clone);
}
