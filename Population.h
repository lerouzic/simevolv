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

#include <iostream>
#include <vector>



class Individual; // This is clearly a design bug. Population.h is called in Individual.h through another way. 


class Population
{
	friend class Canalization;
	friend class Heritability;
	
	public :
	    //constructors
	    Population();
	    Population(long int); // population size
	    Population(const Population&); // copy constructor
	    Population(const std::vector<Individual>&); 
	    Population(const ParameterSet&);
	
	    //operator overload
	    Population& operator = (const Population&); 
	
	    Population reproduce(long int offspr_number = 0) const;
	    
	    void update(void); // ideally, should not be public
	    double mean_phenotype() const; // this is a relic from unidimensional phenotypes. probably used by the fitness routine only
	    
	    long int size() const;

		// Is it safe to leave these two as public? Only friends might be allowed to do this. 
	    void draw_mutation();
	    void make_mutation();
	    	
	    //output
	    void write(int) const;
	    void write_debug(std::ostream&) const;
	    void write_xml(std::ostream&) const;
	    void write_simple(std::ostream&) const;
	    void write_summary(std::ostream&, int) const;
	
	protected :
		// internal functions
	    void initialize(const ParameterSet &);
	    
	    // Stuff for selection
	    std::vector<double> cumul_fitness() const;
	    const Individual & pick_parent(const std::vector<double>&) const;
	    long int search_fit_table(double, const std::vector<double>&) const;
	    long int sequential_search_fit_table(double, const std::vector<double>&) const;
	    long int stl_search_fit_table(double, const std::vector<double>&) const;		
	
	    std::vector<Individual> pop;
	    unsigned int nb_canal_test;
	    unsigned int nb_herit_test;

};

#endif // POPULATION_H_INCLUDED
