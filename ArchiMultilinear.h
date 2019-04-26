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



#ifndef ARCHIMULTILINEAR_H_INCLUDED
#define ARCHIMULTILINEAR_H_INCLUDED

#include "Architecture.h"

#include <iostream>

class Epsilon2
{ 
	public:
		Epsilon2() : dim(0), fixed_dim(false) { }
		Epsilon2(std::size_t d) : dim(d), fixed_dim(true) { }
		Epsilon2(const Epsilon2& c) : dim(c.dim), fixed_dim(c.fixed_dim), value(c.value) { }
		~Epsilon2() { }
		
		void set(std::size_t, std::size_t, allele_type);
		allele_type get(std::size_t, std::size_t) const;
		bool is_init() const;
		
	protected:
		std::size_t dim;
		bool fixed_dim;
		std::vector<allele_type> value;
		
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
			ar & value;
			ar & dim;
			ar & fixed_dim;
		}
        #endif
};

class Epsilon3
{ 
	public:
		Epsilon3() : dim(0), fixed_dim(false) { }
		Epsilon3(std::size_t d) : dim(d), fixed_dim(true) { }
		Epsilon3(const Epsilon3& c) : dim(c.dim), fixed_dim(c.fixed_dim), value(c.value) { }
		~Epsilon3() { }
		
		void set(std::size_t, std::size_t, std::size_t, allele_type);
		allele_type get(std::size_t, std::size_t, std::size_t) const;
		bool is_init() const;
		
	protected:
		std::size_t dim;
		bool fixed_dim;
		std::vector<allele_type> value;
		
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive & ar, const unsigned int version) {
			ar & value;
			ar & dim;
			ar & fixed_dim;
		}
        #endif
};

class ArchiMultilinear : public Architecture
{	
	public :
	    //constructors/destructor
	    ArchiMultilinear(const Architecture&) = delete;
	    ArchiMultilinear(const ParameterSet&);
	    virtual ~ArchiMultilinear() { }
	
	    // getters
	    const Epsilon2& get_epsilon2(std::size_t, std::size_t, std::size_t) const;
	    //Epsilon3& get_epsilon3(std::size_t, std::size_t, std::size_t, std::size_t) const;
	
	    virtual unsigned int nb_phen() const { return nphen; }
	
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo &, bool sdinittest = false, bool sddynamtest = false) const;
	    // virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0, bool test = false) const; // the default is OK

	protected :
	    unsigned int nphen;
	    std::vector<Epsilon2> epsilon2;
	    //std::vector<Epsilon3> epsilon3;
	    
	protected:
		ArchiMultilinear() { }; // Necessary for serialization
		
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
			ar & boost::serialization::base_object<Architecture>(*this);
			ar & nphen;
			ar & epsilon2;
			//ar & epsilon3;
		}
        #endif
};

#endif // ARCHIMULTILINEAR_H_INCLUDED
