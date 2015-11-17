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



#ifndef ARCHITECTURE_H_INCLUDED
#define ARCHITECTURE_H_INCLUDED

#include "Parameters.h"
#include "GeneticMap.h"
#include "Allele.h"
#include "Genotype.h"
#include "Phenotype.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class Architecture  	/* Pure virtual class */
{
	public :
	    //constructors/destructor
	    Architecture(const Architecture&) = delete;
	    Architecture (const ParameterSet&); 
	    
	    virtual ~Architecture();
	
	    // instance / initialization
	    static void initialize(const ParameterSet&); // from the parameter file
	    static void initialize(const std::string&);  // from the architecture file
	    
	    static /* const */ Architecture* Get();
	    static void terminate();
	
	    // getters
	    unsigned int nb_loc() const;
	    unsigned int all_size() const;
	    double mutation_rate(unsigned int) const;
	    double mutation_sd(unsigned int) const;
	    double mutation_sd_test(unsigned int) const;
	    double recombination_rate(unsigned int) const;
	    		
		// to be defined by inherited classes 
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir) const = 0; // no default
	    virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int loc = 0) const;
	    virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
	    virtual std::shared_ptr<Allele> allele_mutation_test(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
	
	protected :
	    static Architecture* instance;

	    GeneticMap gmap;
	    unsigned int nloc; // number of loci
	    unsigned int sall; // size of alleles
	    std::vector<double> mutrate;
	    std::vector<double> mutsd;
	    std::vector<double> mutsd_test;
	    
	    std::string iofile; // this should not be serialized
	    
		Architecture() {} // default constructor necessary for serialization
	
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);
};

// Templates have to be written in the header file
template<class Archive>
void Architecture::serialize(Archive & ar, const unsigned int version) 
{
	ar & gmap;
	ar & nloc;
	ar & sall; 
	ar & mutrate;
	ar & mutsd;
	ar & mutsd_test;	
}

#endif // ARCHITECTURE_H_INCLUDED
