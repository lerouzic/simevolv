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
#include "Haplotype.h"
#include "Genotype.h"
#include "Phenotype.h"

#include <iostream>
#include <vector>
#include <memory>



class Architecture  	/* Pure virtual class */
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
	
	    //functions
	    unsigned int nb_loc() const;
	    unsigned int all_size() const;
	    double mutation_rate(unsigned int) const;
	    double mutation_sd(unsigned int) const;
	    double recombination_rate(unsigned int) const;
	    		
		//inheritance
	    virtual Phenotype phenotypic_value(const Genotype&) const = 0;
	    virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int loc = 0) const;
	    virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0) const;

	
	protected :
	    GeneticMap gmap;
	    unsigned int nloc;
	    unsigned int sall;
	    std::vector<double> mutrate;
	    std::vector<double> mutsd;

};


#endif // ARCHITECTURE_H_INCLUDED

