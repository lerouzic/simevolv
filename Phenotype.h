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



#ifndef PHENOTYPE_H_INCLUDED
#define PHENOTYPE_H_INCLUDED

#include "Statistics.h"

#include <iostream>
#include <vector>



class Phenotype
{
	public:
		//constructors/destructor
		Phenotype();
		Phenotype(const double);
		Phenotype(const std::vector<double> &);
		Phenotype(const Phenotype &);
		~Phenotype();
		
		//initialization
		void initialize(const std::vector<double>&);
		
		//operator overload
		Phenotype& operator= (const Phenotype &);
		double operator[] (const unsigned int index) const;
				
		//functions
		void copy(const Phenotype &);
		unsigned int dimensionality() const;
		
		//output
		void write_debug (std::ostream&) const;	
		void write_simple (std::ostream&) const;	
	    friend std::ostream& operator << (std::ostream&, const Phenotype&); // probably doesn't need to be friend
		
	protected:
		std::vector<double> pheno;
};




///////////// Specific class to get multivariate statistics on phenotypes. 

class PhenotypeStat: public MultivariateStat
{
	public:
		//constructor
		PhenotypeStat(const std::vector<Phenotype> &);
		
		//functions
		Phenotype means_phen() const {return(Phenotype(means()));}
		Phenotype vars_phen() const {return(Phenotype(vars()));}
		unsigned int dimensionality() const {return(data.size());}
		static std::vector<std::vector<double> > transpose_phen_matrix(const std::vector<Phenotype> &);
		
	protected:
		
};


#endif
