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



#ifndef GENOTYPE_H_INCLUDED
#define GENOTYPE_H_INCLUDED

#include "types.h"
#include "Haplotype.h"
#include "Parameters.h"
#include "Phenotype.h"

#include <iostream>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#endif

class Genotype
{
	friend class Population;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiMasel;
    friend class ArchiSiegal; 
    friend class ArchiM2;
    friend class ArchiBoolean;

	public:
		virtual ~Genotype() {}
		virtual Genotype* clone() const = 0;
	
	    unsigned int nb_loc() const;
	    unsigned int all_size() const;
	    virtual unsigned int ploidy() const = 0;
	    virtual void draw_mutation() = 0;
	    virtual void make_mutation(bool test = false) = 0;
		virtual Haplotype produce_gamete() const = 0;
		
		virtual std::vector<allele_type> combine_at_loc(unsigned int, std::vector<allele_type> (*combineFUN)(const Allele &, const Allele &)) const = 0;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            // nothing to do
        }
        #endif
};

class DiploGenotype : public Genotype {
	public:
        DiploGenotype() { }
		DiploGenotype(const DiploGenotype &);
		DiploGenotype(const Haplotype & gam_father, const Haplotype & gam_mother);
		DiploGenotype(const ParameterSet &);
		DiploGenotype* clone() const;

		unsigned int ploidy() const {return 2;}
	    void draw_mutation();
	    void make_mutation(bool test = false);
		Haplotype produce_gamete() const;
		
		std::vector<allele_type> combine_at_loc(unsigned int, std::vector<allele_type> (*combineFUN)(const Allele &, const Allele &)) const;	
		
	protected:
	    Haplotype gam_father;
	    Haplotype gam_mother;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & boost::serialization::base_object<Genotype>(*this);
            ar & gam_father;
            ar & gam_mother;
        }        
        #endif
};

class HaploGenotype : public Genotype {
	public:
        HaploGenotype() { }
		HaploGenotype(const HaploGenotype &);
		HaploGenotype(const Haplotype & gam_father, const Haplotype & gam_mother);
		HaploGenotype(const ParameterSet &);
		HaploGenotype* clone() const;
		
		unsigned int ploidy() const {return 1;}
	    void draw_mutation();
	    void make_mutation(bool test = false);
		Haplotype produce_gamete() const;		
		
		std::vector<allele_type> combine_at_loc(unsigned int, std::vector<allele_type> (*combineFUN)(const Allele &, const Allele &)) const;
			
	protected: 
		Haplotype gam;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & boost::serialization::base_object<Genotype>(*this);
            ar & gam;
        }           
        #endif
};

#endif // GENOTYPE_H_INCLUDED
