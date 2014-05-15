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
		Phenotype(const std::vector<double> &, const std::vector<double> &);
		Phenotype(const Phenotype &);
		~Phenotype();
		
		//initialization
		void initialize(const std::vector<double>&);
		
		//operator overload
		Phenotype& operator= (const Phenotype &);
		
		// getters
		double operator[] (const unsigned int index) const;
		double get_pheno(const unsigned int index) const;
		double get_unstab(const unsigned int index) const;
		void add_pheno(const unsigned int index, const double effect);
				
		// number of phenotypes
		unsigned int dimensionality() const;
		
		//output
		void write_debug (std::ostream&) const;	
		void write_simple (std::ostream&) const;	
	    friend std::ostream& operator << (std::ostream&, const Phenotype&); 
		
	protected:
		std::vector<double> pheno;
		std::vector<double> unstabpheno;
};




///////////// Specific class to get multivariate statistics on phenotypes. 

class PhenotypeStat
{
	public:
		//constructor
		PhenotypeStat(const std::vector<Phenotype> &);
		
		//functions
		std::vector<double> means_phen() const {return(pheno.means());}
		std::vector<double> vars_phen() const {return(pheno.vars());}
		
		std::vector<double> means_unstab() const {return(unstab.means());}
		std::vector<double> vars_unstab() const {return(unstab.vars());}
		
		unsigned int dimensionality() const {return(pheno.means().size());}
		static std::vector<std::vector<double> > transpose_phen_matrix(const std::vector<Phenotype> &);
		static std::vector<std::vector<double> > transpose_unstabphen_matrix(const std::vector<Phenotype> &);
		
	protected:
		MultivariateStat pheno;
		MultivariateStat unstab;
		
};

class InvertedMStat: public MultivariateStat
{
	public:
		InvertedMStat(const std::vector<std::vector<double> > &);
	
	protected:
	// functions
	static std::vector<std::vector<double> > transpose_double_matrix(const std::vector<std::vector<double> > &);
};


#endif
