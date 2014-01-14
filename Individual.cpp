// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Individual.h"
#include "main.h"
#include "Fitness.h"
#include "Environment.h"
#include "Architecture.h"

#include <vector>

using namespace std;



// constructors and destructor
Individual::Individual()
    : genotype()
{
    initialize();
}


Individual::Individual(const Haplotype& gam_father, const Haplotype& gam_mother)
    : genotype(gam_father, gam_mother)
{
    initialize();
}


Individual::Individual(const Individual& copy)
    : genotype(copy.genotype)
    , genot_value(copy.genot_value)
    , phenotype(copy.phenotype)
    , fitness(0.0)
{
}


Individual::Individual(const ParameterSet& param)
{
    initialize();
}


Individual::~Individual()
{
}


// operator overload

Individual & Individual::operator= (const Individual& copy)
{
    if (this == &copy)
        return (*this);

    genotype=copy.genotype;
    genot_value=copy.genot_value;
    phenotype=copy.phenotype;
    fitness=copy.fitness;

    return(*this);
}


// instance and initialization

void Individual::initialize()
{
    Architecture * archi = Architecture::Get();
    genot_value = archi -> phenotypic_value(genotype);
    phenotype = Environment::rand_effect(genot_value);
    fitness = 0;
}


// functions

void Individual::update_fitness(const Population & pop)
{
    fitness = Fitness::compute(phenotype, pop);
}


void Individual::update_fitness(const double population_value)
{
    fitness = Fitness::compute(phenotype, population_value);
}


double Individual::get_fitness() const
{
    return(fitness);
}

Phenotype Individual::get_genot_value() const
{
	return(genot_value);
}

Phenotype Individual::get_phenotype() const
{
    return(phenotype);
}


Individual Individual::mate(const Individual& father, const Individual& mother)
{
    Individual offspring(father.produce_gamete(), mother.produce_gamete());
    return(offspring);
}


Haplotype Individual::produce_gamete() const
{
    Haplotype gamete(genotype.recombine());
    gamete.draw_mutation();
    return(gamete);
}


void Individual::draw_mutation()
{
    genotype.draw_mutation();
    genot_value = Architecture::Get() -> phenotypic_value(genotype);
}


void Individual::make_mutation()
{
    genotype.make_mutation();
    genot_value = Architecture::Get() -> phenotypic_value(genotype);
    // phenotype = Environment::rand_effect(genot_value); // This step is not obvious
    // fitness = 0; // computing fitness requires additional information, call update_fitness with the proper argument
}

Individual Individual::test_canalization(unsigned int nb_mut, const Population & pop) const 
{
	Individual clone(*this);
	for (unsigned int mut = 0; mut < nb_mut; mut++) {
		clone.make_mutation();
	}
	clone.update_fitness(pop);
	return(clone);		
}

// output

void Individual::write_debug(ostream& out) const
{
    out << "Genotype:" << endl;
    genotype.write_debug(out);
    phenotype.write_debug(out);
    out << "Fitness:\t" << fitness << endl;
}


void Individual::write_xml(ostream& out) const
{
    out << "xml output: not implemented yet.\n";
}


void Individual::write_simple(ostream& out) const
{
    out << "Pheno:";
    phenotype.write_simple(out);
    out << "\tFitness: " << get_fitness() << endl;
}
