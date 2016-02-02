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



#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include "Phenotype.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "Parameters.h"
#include "EpigeneticInfo.h"

#include <iostream>
#include <string>


class Population;

class Individual
{
	friend class Population;	
	public :
	    // constructors/destructor
	    Individual() = delete;
	    Individual(const Haplotype&, const Haplotype&);
	    Individual(const Haplotype&, const Haplotype&, const EpigeneticInfo&);
	    Individual(const Individual&);
	    Individual(const ParameterSet&);
	    virtual ~Individual();
	
	    // operator overload
	    Individual & operator = (const Individual&);
	    int operator == (const Individual&) const;
	
	    // instance/initialization
	    void initialize();
	    void update_fitness(const Population &);

		// getters
	    double get_fitness() const;
	    double get_epigenet() const;
	    Phenotype get_phenotype() const;
	    std::string write_debug(unsigned int) const;
	    EpigeneticInfo make_epiinfo() const;
	    
	    // reproduction
	    Haplotype produce_gamete() const;
	    static Individual mate(const Individual&, const Individual&);
	    
	    void draw_mutation();
	    void make_mutation(bool test = false);
	    	    		
	    /* Make a clone with some mutations
		- the first parameter is the number of mutations */	    	    		
	    Individual test_canalization(unsigned int, const Population &) const;  
	
	protected :
	    Genotype genotype;
	    Phenotype phenotype; // phenotype + environmental effect 
	    double fitness;
	    double epigenet; // to be transmitted to the offspring
	    
	    EpigeneticInfo epiinfo; // from the mother
};

#endif // INDIVIDUAL_H_INCLUDED
