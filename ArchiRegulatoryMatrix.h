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
#include <boost/numeric/ublas/vector.hpp>

#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>

class ArchiRegulatoryMatrix : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatoryMatrix(const ArchiRegulatoryMatrix&) = delete;  
	    ArchiRegulatoryMatrix(const ParameterSet&);
	    virtual ~ArchiRegulatoryMatrix() = 0; // So that the class is abstract
		
		// functions 
		virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
		virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo&, bool sdinittest = false, bool sddynamtest = false) const;
	
	protected :
		unsigned int sall;
		std::vector<double> So;
		double recur;
		std::vector<std::vector<double>> connectivity_matrix; // this contains initial allelic values (for clonal pops), not only 0 or 1
		unsigned int timesteps;
		unsigned int calcsteps;
				
	    //functions
		virtual double sigma(double h) const;
		void init_connectivity_matrix(const ParameterSet &);		
		virtual void haircut(boost::numeric::ublas::vector<double> & vec) const;
		virtual void haircut(std::vector<double> & vec) const;
		virtual void haircut(double &) const;
		
	protected:
		ArchiRegulatoryMatrix() {}
		
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);			
};

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

class ArchiWagner : public ArchiRegulatoryMatrix
{
	public :
	    //constructors/destructor
	    ArchiWagner(const ArchiWagner &) = delete;
	    ArchiWagner(const ParameterSet&);
	    virtual ~ArchiWagner();
		
	protected :
		// Inherited functions
		virtual double sigma(double h) const;
		virtual void haircut(double &) const;
		
	protected:
		ArchiWagner() {}
		
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);					
};

template<class Archive>
void ArchiWagner::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
}

class ArchiSiegal : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
		ArchiSiegal(const ArchiSiegal &) = delete;
	    ArchiSiegal(const ParameterSet&);
	    virtual ~ArchiSiegal();

	protected :
		double basal;
		
		// Inherited functions
		virtual double sigma(double h) const; 
		virtual void haircut(double &) const;	
		
	protected:
		ArchiSiegal() {}
		
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);					
};

template<class Archive>
void ArchiSiegal::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
	ar & basal;
}		


class ArchiM2 : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiM2(const ArchiM2 &) = delete;
	    ArchiM2(const ParameterSet&);
	    virtual ~ArchiM2();
			
	protected :
		double basal;
		
		// Inherited functions
		virtual double sigma(double h) const; 
		virtual void haircut(double &) const;
		
	protected:
		ArchiM2() {}

	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);			
};

template<class Archive>
void ArchiM2::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<ArchiRegulatoryMatrix>(*this);	
	ar & basal;
}

#endif // ARCHIREGULATORYMATRIX_H_INCLUDED
