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



#ifndef ARCHIREGULATORYMATRIX_H_INCLUDED
#define ARCHIREGULATORYMATRIX_H_INCLUDED

#include "Architecture.h"

#include <cmath>
#include <iostream>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#endif

class ArchiRegulatoryMatrix : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatoryMatrix(const ArchiRegulatoryMatrix&) = delete;  
	    ArchiRegulatoryMatrix(const ParameterSet&);
	    virtual ~ArchiRegulatoryMatrix() = 0; // So that the class is abstract
		
		// functions 
        virtual unsigned int nb_phen() const;
        
		virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
	    virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0, bool test = false) const;

		virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo&, bool sdinittest = false, bool sddynamtest = false) const;
	
	protected :
		std::vector<pheno_type> So;
		rate_type recur;
		std::vector<std::vector<allele_type>> connectivity_matrix; // this contains initial allelic values (for clonal pops), not only 0 or 1
		unsigned int timesteps;
		unsigned int calcsteps;
				
	    //functions
        virtual void sigma(pheno_type &) const;
        virtual void haircut(pheno_type &) const;
        
        // vectorized versions for optimization
        virtual void sigma_v(std::vector<pheno_type>&) const; 
		virtual void haircut_v(std::vector<pheno_type>&) const;
        virtual void plasticity_v(std::vector<pheno_type>&) const;
        virtual void recur_v(std::vector<pheno_type>&, const std::vector<pheno_type>&) const;
        virtual void enviro_v(std::vector<pheno_type>&, bool) const;
        
        
		void init_connectivity_matrix(const ParameterSet &);		
		
	protected:
		ArchiRegulatoryMatrix() {}
		
	private:
    #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);			
    #endif
};

#ifdef SERIALIZATION_TEXT
template<class Archive>
void ArchiRegulatoryMatrix::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<Architecture>(*this);	
	ar & sall;
	ar & So;
	ar & recur;
	ar & connectivity_matrix;
	ar & timesteps;
	ar & calcsteps;
}
#endif

class ArchiWagner : public ArchiRegulatoryMatrix
{
	public :
	    //constructors/destructor
	    ArchiWagner(const ArchiWagner &) = delete;
	    ArchiWagner(const ParameterSet&);
	    virtual ~ArchiWagner();
		
	protected :
		// Inherited functions
		virtual void sigma(pheno_type&) const;
        virtual void sigma_v(std::vector<pheno_type>&) const;
		virtual void haircut(pheno_type &) const;
        virtual void haircut_v(std::vector<pheno_type>&) const;
		
	protected:
		ArchiWagner() {}
		
	private:
    #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);					
    #endif
};

#ifdef SERIALIZATION_TEXT
template<class Archive>
void ArchiWagner::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
}
#endif

class ArchiSiegal : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
		ArchiSiegal(const ArchiSiegal &) = delete;
	    ArchiSiegal(const ParameterSet&);
	    virtual ~ArchiSiegal();

	protected :
		pheno_type basal;
		
		// Inherited functions
		virtual void sigma(pheno_type&) const; 
        virtual void sigma_v(std::vector<pheno_type>&) const;
		virtual void haircut(pheno_type&) const;	
        virtual void haircut_v(std::vector<pheno_type>&) const;
		
	protected:
		ArchiSiegal() {}
		
	private:
    #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);					
    #endif
};

#ifdef SERIALIZATION_TEXT
template<class Archive>
void ArchiSiegal::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
	ar & basal;
}		
#endif

class ArchiM2 : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiM2(const ArchiM2 &) = delete;
	    ArchiM2(const ParameterSet&);
	    virtual ~ArchiM2();
			
	protected :
		pheno_type basal;
        pheno_type b1; // = 1/basal - 1
        pheno_type b2; // = 1/(basal*(basal-1))
		
		// Inherited functions
		virtual void sigma(pheno_type&) const; 
        virtual void sigma_v(std::vector<pheno_type>&) const;
		virtual void haircut(pheno_type&) const;
        virtual void haircut_v(std::vector<pheno_type>&) const;
		
	protected:
		ArchiM2() { };

	private:
    #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);			
    #endif
};

#ifdef SERIALIZATION_TEXT
template<class Archive>
void ArchiM2::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
	ar & basal;
    ar & b1;
    ar & b2;
}
#endif

#endif // ARCHIREGULATORYMATRIX_H_INCLUDED
