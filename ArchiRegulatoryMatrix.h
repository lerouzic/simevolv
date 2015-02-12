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

#include <math.h>
#include <iostream>


class ArchiRegulatoryMatrix : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatoryMatrix() = delete; 	// should not be used (C++11)
	    ArchiRegulatoryMatrix(const ArchiRegulatoryMatrix&) = delete;  
	    ArchiRegulatoryMatrix(const ParameterSet&);
	    virtual ~ArchiRegulatoryMatrix() {}
		
		// functions 
		virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
		virtual Phenotype phenotypic_value(const Genotype&) const;
	
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
};

class ArchiWagner : public ArchiRegulatoryMatrix
{
	public :
	    //constructors/destructor
	    ArchiWagner() = delete;
	    ArchiWagner(const ArchiWagner &) = delete;
	    ArchiWagner(const ParameterSet&);
	    virtual ~ArchiWagner() {}
		
	protected :
		// Inherited functions
		virtual double sigma(double h) const;
};	

class ArchiSiegal : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiSiegal() = delete;
		ArchiSiegal(const ArchiSiegal &) = delete;
	    ArchiSiegal(const ParameterSet&);
	    virtual ~ArchiSiegal() {}

	protected :
		double basal;
		
		// Inherited functions
		virtual double sigma(double h) const; 
};

class ArchiM2 : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiM2() = delete;
	    ArchiM2(const ArchiM2 &) = delete;
	    ArchiM2(const ParameterSet&);
	    virtual ~ArchiM2() {}
			
	protected :
		double basal;
		
		// Inherited functions
		virtual double sigma(double h) const; 
};

#endif // ARCHIREGULATORYMATRIX_H_INCLUDED
