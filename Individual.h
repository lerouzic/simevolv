// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
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

#include "Population.h"
#include "Phenotype.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "Parameters.h"

#include <iostream>



class Population;


class Individual
{
	public :
	    // constructors/destructor
	    Individual(const Haplotype&, const Haplotype&);
	    Individual(const Individual&);
	    Individual(const ParameterSet&);
	    virtual ~Individual();
	
	    // operator overload
	    Individual & operator = (const Individual&);
	    int operator == (const Individual&) const;
	
	    // instance/initialization
	    void initialize();
	
	    // functions
	    void update_fitness(const Population &);
	    double get_fitness() const;
	    Phenotype get_phenotype() const;
	    Phenotype get_genot_value() const;
	    
	    Haplotype produce_gamete() const;
	    static Individual mate(const Individual&, const Individual&);
	    void draw_mutation();
	    void make_mutation();
	    	    		
	    Individual test_canalization(unsigned int, const Population &) const;  /* the first parameter is the number of mutations */
	
	    //output
	    void write_debug (std::ostream&) const;
	    void write_xml (std::ostream&) const;
	    void write_simple(std::ostream&) const;
	
	
	public :
	    Genotype genotype;
	    Phenotype genot_value; // a bit counterintutive, and problematic for the future. 
	    Phenotype phenotype;
	    double fitness;
};

#endif // INDIVIDUAL_H_INCLUDED
