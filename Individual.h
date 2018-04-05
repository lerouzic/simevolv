// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2017 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>
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

#include "types.h"
#include "Phenotype.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "Parameters.h"
#include "EpigeneticInfo.h"
#include "Fitness.h"

#include <iostream>
#include <string>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unique_ptr.hpp>
#endif

class Population;

class Individual
{
	friend class Population;	
	public :
	    // constructors/destructor
	    Individual();
	    Individual(const Individual&);
	    Individual(const ParameterSet&);
                
	    Individual(const Haplotype&, const Haplotype&, const unsigned int);
	    Individual(const Haplotype&, const Haplotype&, const unsigned int, const EpigeneticInfo&);
        Individual(const Genotype&);
        Individual(const Genotype&, const EpigeneticInfo&);

	    virtual ~Individual();
	
	    // operator overload
	    Individual & operator = (const Individual&);
	    int operator == (const Individual&) const;
	
	    // instance/initialization        
	    void update_fitness(const Population &);
        void update_phenotype();
        void update_epigenet(rate_type);
	    // void initialize(); // Updates whatever possible (useless?)

		// getters
	    fitness_type get_fitness() const;
	    rate_type get_epigenet() const;
	    Phenotype get_phenotype() const;
	    unsigned int ploidy() const;
	    
	    EpigeneticInfo make_epiinfo() const;
	    
	    // reproduction
	    Haplotype produce_gamete() const;
	    static Individual mate(const Individual&, const Individual&);
        Individual clone() const;
	    
	    void draw_mutation();
	    void make_mutation(bool test = false);
	    	    		
	    /* Make a clone with some mutations
		- the first parameter is the number of mutations */	    	    		
	    Individual test_canalization(unsigned int, const Population &) const;
	    
	    Individual test_disturb(const Population &) const;
	    
	    Individual test_enviro(const Population &) const;
	
	protected :
	    std::unique_ptr<Genotype> genotype;
	    Phenotype phenotype; 
	    fitness_type fitness;
	    rate_type epigenet; // to be transmitted to the offspring
	    
	    EpigeneticInfo epiinfo; // from the mother
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & genotype;
            ar & phenotype;
            ar & fitness;
            ar & epigenet;
            ar & epiinfo;            
        } 
        #endif
};

#endif // INDIVIDUAL_H_INCLUDED
