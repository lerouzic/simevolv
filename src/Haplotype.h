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



#ifndef HAPLOTYPE_H_INCLUDED
#define HAPLOTYPE_H_INCLUDED

#include "Allele.h"
#include "Parameters.h"

#include <iostream>
#include <string>
#include <vector>
#include <memory> // shared_ptr

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/shared_ptr.hpp>
#endif

class Haplotype
{
    friend class Individual;
	friend class Population;
    friend class Genotype;
    friend class HaploGenotype;
    friend class DiploGenotype;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiFKL;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiSiegal; 
    friend class ArchiM2;
    friend class ArchiBoolean;

	public :
	    //constructors/destructor
        Haplotype() { }
	    Haplotype(const ParameterSet &);
	    Haplotype(const Haplotype &);
	    Haplotype(const std::vector<std::shared_ptr<Allele> > &);
		~Haplotype() {}
	
	    //operator overload
	    int operator == (const Haplotype&) const;
	    int operator != (const Haplotype&) const;
	
	    //functions
	    unsigned int nb_loc() const;
	    void draw_mutation();				
	    void make_mutation(bool test = false); 				
	    void make_mutation(unsigned int, bool test = false);
	    void make_mutmutation();
	    void make_mutmutation(unsigned int);
	    
		std::shared_ptr<const Allele> allele_at_loc(unsigned int) const;
	    
	    static Haplotype recombine(const Haplotype &, const Haplotype &);
	    
	    std::string write_debug() const;
	
	protected :
	    std::vector<std::shared_ptr<Allele> > haplotype;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & haplotype;
        }
        #endif        
};

#endif // HAPLOTYPE_H_INCLUDED
