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

#include "types.h"
#include "Parameters.h"
#include "Phenotype.h"
#include "Individual.h"
#include "Fitness.h"

#include <iostream>
#include <vector>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

class Population
{
	friend class GeneticCanalization;
	friend class DisturbCanalization;
	friend class EnviroCanalization;
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
	
		void update_param(const ParameterSet&);
	    Population reproduce(long int offspr_number = 0) const;
	    
	    Phenotype mean_phenotype() const; // Probably used by the Fitness class only
	    long int size() const;

	    void draw_mutation();
	    void make_mutation();
	    	
	    // input / output
            // write summary
	    void write(std::ostream &, int) const;
            // serialization
        #ifdef SERIALIZATION_TEXT
        friend std::ostream& operator << (std::ostream&, const Population&);
        friend std::istream& operator >> (std::istream&, Population&);
        #endif
        
	protected :
		// internal functions
	    void initialize(const ParameterSet &);
	    void update(void); 
	    
	    // Stuff for selection
	    std::vector<fitness_type> cumul_fitness() const;
	    const Individual & pick_parent(const std::vector<fitness_type>&) const;
	    // different algorithms to optimize weighted random picking of parents
	    long int search_fit_table(fitness_type, const std::vector<fitness_type>&) const;
	    long int sequential_search_fit_table(fitness_type, const std::vector<fitness_type>&) const;
	    long int stl_search_fit_table(fitness_type, const std::vector<fitness_type>&) const;		
	
	    std::vector<Individual> pop;
	    
	    // output parameters, to be stored and copied (design problem?)
	    rate_type selfing_rate;	 
        rate_type clonal_rate;   
	    unsigned int nb_canal_test;
	    unsigned int nb_herit_test;
	    unsigned int nb_direpi_test;
	    std::string out_geno;
	    std::string out_unstab;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version){
            ar & pop;
            ar & selfing_rate;
            ar & clonal_rate;
            ar & nb_canal_test;
            ar & nb_herit_test;
            ar & nb_direpi_test;
            ar & out_geno;
            ar & out_unstab;
        }     
        #endif
};


#endif // POPULATION_H_INCLUDED
