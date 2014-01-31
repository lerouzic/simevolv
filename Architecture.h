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



#ifndef ARCHITECTURE_H_INCLUDED
#define ARCHITECTURE_H_INCLUDED

#include "Parameters.h"
#include "GeneticMap.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"

#include <iostream>
#include <vector>



class Architecture
{
	public :
	    //constructors/destructor
	    Architecture();
	    Architecture(const Architecture&);
	    Architecture (const ParameterSet&);
	    virtual ~Architecture(){};
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    // instance / initialization
	    static Architecture* instance;
	    static void initialize(const ParameterSet&);
	    static Architecture* Get();
	    static Architecture* Get(const ParameterSet*);
	    static Architecture* Get(const ParameterSet&);
	
	    //functions
	    int nb_loc() const;
	    int all_size() const;
	    double init_all() const;
	    double mutation_rate(int) const;
	    double mutation_sd(int) const;
	    double recombination_rate(int) const;
	    void draw_mutation(const Haplotype&) const;
	    double make_mutation(int, std::vector<Allele>) const;
		
		//inheritance
	    virtual double phenotypic_value(const Genotype&) const;
	
	protected :
	    GeneticMap gmap;
	    int nloc;
	    int sall;
	    double iall;
	    std::vector<double> mutrate;
	    std::vector<double> mutsd;

};


#endif // ARCHITECTURE_H_INCLUDED

