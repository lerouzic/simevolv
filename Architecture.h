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

#include "types.h"
#include "Parameters.h"
#include "GeneticMap.h"
#include "Allele.h"
#include "Genotype.h"
#include "Phenotype.h"
#include "EpigeneticInfo.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

class Architecture  	/* Pure virtual class */
{
	public :
	    //constructors/destructor
	    Architecture(const Architecture&) = delete;
	    Architecture (const ParameterSet&); 
	    
	    virtual ~Architecture();
	
	    // instance / initialization
	    static void initialize(const ParameterSet&); // from the parameter file
	    static void update_param(const ParameterSet&);
	    
	    static void load(const std::string&);  // from the architecture file
        static void save(const std::string&);         
        
	    static /* const */ Architecture* Get();	
  
	    // getters
	    unsigned int nb_loc() const;
        virtual unsigned int nb_phen() const;
	    unsigned int all_size() const;
	    rate_type mutation_rate(unsigned int) const;
	    allele_type mutation_sd(unsigned int) const;
	    allele_type mutation_sd_test(unsigned int) const;
	    rate_type recombination_rate(unsigned int) const;
	    		
		// to be defined by inherited classes 
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo&, bool sdinittest = false, bool sddynamtest = false) const = 0; // no default
	    virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int loc = 0) const;
	    virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
	    virtual std::shared_ptr<Allele> allele_mutation_test(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
	
	protected :
	    static Architecture* instance;

	    GeneticMap gmap;
	    unsigned int nloc; // number of loci
	    unsigned int sall; // size of alleles
	    std::vector<rate_type> mutrate;
	    std::vector<allele_type> mutsd;
	    std::vector<allele_type> mutsd_test;
        
        std::vector<pheno_type> plasticity_strength;
        std::vector<pheno_type> plasticity_signal;
	    	    
		Architecture() {} // default constructor necessary for serialization
	    virtual void update_param_internal(const ParameterSet&);		
	
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
        {
            ar & gmap;
            ar & nloc;
            ar & sall; 
            ar & mutrate;
            ar & mutsd;
            ar & mutsd_test;
            ar & plasticity_strength;
            ar & plasticity_signal;
        }
        #endif
};

#endif // ARCHITECTURE_H_INCLUDED
