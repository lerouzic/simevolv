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



#ifndef POPULATION_H_INCLUDED
#define POPULATION_H_INCLUDED

#include "Parameters.h"
#include "Phenotype.h"
#include "Individual.h"

#include <iostream>
#include <vector>


class Population
{
	friend class Canalization;
	friend class Heritability;
	friend class Direpistasis;

	public :
	    //constructors
	    Population();
	    Population(const Population&); // copy constructor
	    Population(const std::vector<Individual>&); 
	    Population(const ParameterSet&);
	
	    //operator overload
	    Population& operator = (const Population&); 
	
	    Population reproduce(long int offspr_number = 0) const;
	    
	    Phenovec mean_phenotype() const; // Probably used by the Fitness class only
	    long int size() const;

	    void draw_mutation();
	    void make_mutation();
	    	
	    //output
	    void write(std::ostream &, int) const;
	    void write_debug(std::ostream &) const;
	
	protected :
		// internal functions
	    void initialize(const ParameterSet &);
	    void update(void); 
	    
	    // Stuff for selection
	    std::vector<double> cumul_fitness() const;
	    const Individual & pick_parent(const std::vector<double>&) const;
	    // different algorithms to optimize weighted random picking of parents
	    long int search_fit_table(double, const std::vector<double>&) const;
	    long int sequential_search_fit_table(double, const std::vector<double>&) const;
	    long int stl_search_fit_table(double, const std::vector<double>&) const;		
	
	    std::vector<Individual> pop;
	    unsigned int nb_canal_test;
	    unsigned int nb_herit_test;
	    unsigned int nb_direpi_test;
	    std::string out_geno;
	    std::string out_unstab;
};

#endif // POPULATION_H_INCLUDED
