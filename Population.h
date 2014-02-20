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



#ifndef POPULATION_H_INCLUDED
#define POPULATION_H_INCLUDED

#include "Parameters.h"
#include "Individual.h"
#include "Canalization.h"

#include <iostream>
#include <vector>


class Individual;

class Population
{
	friend class Canalization;

	public :
	    //constructors/destructors
	    Population();
	    Population(long int);
	    Population(const Population&);
	    Population(const std::vector<Individual>&);
	    Population(const ParameterSet&);
	
	    //operator overload
	    Population& operator = (const Population&);
	    int operator == (const Population&) const;
	
	    //instance/initialization
	    void initialize(const ParameterSet &);
	
	    //functions
	    Population reproduce(long int offspr_number = 0) const;
	    void update(void);
	    //~ std::vector<double> phenotypes() const;
	    double mean_phenotype() const;
	    long int size() const;
	    std::vector<double> cumul_fitness() const;
	    const Individual & pick_parent(const std::vector<double>&) const;
	    long int search_fit_table(double, const std::vector<double>&) const;
	    long int sequential_search_fit_table(double, const std::vector<double>&) const;
	    long int stl_search_fit_table(double, const std::vector<double>&) const;
	    void draw_mutation();
	    void make_mutation();
	    	
	    //output
	    void write(int) const;
	    void write_debug(std::ostream&) const;
	    void write_xml(std::ostream&) const;
	    void write_simple(std::ostream&) const;
	    void write_summary(std::ostream&, int) const;
	
	
	protected :
	    std::vector<Individual> pop;
	    unsigned int nb_canal_test;

};

#endif // POPULATION_H_INCLUDED
