// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
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



class ArchiRegulatoryMatrix : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatoryMatrix();
	    ArchiRegulatoryMatrix(const ParameterSet&);
	    ~ArchiRegulatoryMatrix() {}
		
		// functions 
		virtual double sigma(double h) const {return (h);}
		virtual Phenotype phenotypic_value(const Genotype&) const;
		std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
	
	protected :
		unsigned int sall;
		std::vector<double> So;
		std::vector<std::vector<double> > connectivity_matrix; // this contains initial allelic values (for clonal pops), not only 0 or 1
		unsigned int timesteps;
		unsigned int calcsteps; 
				
	    //functions
		void init_connectivity_matrix(const ParameterSet &);		
};


class ArchiWagner : public ArchiRegulatoryMatrix
{
	public :
	    //constructors/destructor
	    ArchiWagner() {assert(false);}
	    ArchiWagner(const ParameterSet&);
	    ~ArchiWagner() {}
		
		// Inherited functions
		double sigma(double h) const {if (h<0){return (-1.);} else if (h>0) {return (1.);} else {return (0.);}}	
};	



class ArchiMasel : public ArchiRegulatoryMatrix
{
	public :
	    //constructors/destructor
	    ArchiMasel() {assert(false);}
	    ArchiMasel(const ParameterSet&);
	    ~ArchiMasel() {}
		
		// Inherited functions
		double sigma(double h) const {if (h>=0){return (1.);} else {return (0.);}}	
};


class ArchiSiegal : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiSiegal() {assert(false);}
	    ArchiSiegal(const ParameterSet&);
	    ~ArchiSiegal() {}
		
		// Inherited functions
		double sigma(double h) const {return ((2 / (1 + exp(-basal*h)) ) -1);}
		
	protected :
		double basal;
};


class ArchiM2 : public ArchiRegulatoryMatrix
{
	public : 
		//constructors/destructor
	    ArchiM2() {assert(false);}
	    ArchiM2(const ParameterSet&);
	    ~ArchiM2() {}
		
		// Inherited functions
		double sigma(double h) const {return (1 / (1 + exp((-h/(basal*(1-basal)))+log(1/basal-1)) ));}
	
	protected :
		double basal;
};



#endif // ARCHIREGULATORYMATRIX_H_INCLUDED
