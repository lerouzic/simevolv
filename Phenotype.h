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

#include <boost/serialization/serialization.hpp>

class Phenovec;
void outformat(std::ostream &, const Phenovec &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");
class Phenotype;
void outformat(std::ostream &, const Phenotype &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");


// Phenovec is "just" a vector of double, storing phenotypic values for several characters
class Phenovec {
	friend std::ostream& operator << (std::ostream&, const Phenovec &);
	friend void outformat(std::ostream &, const Phenovec &, 
		unsigned int width, unsigned int precision, std::string sep);    
        
    public:
        Phenovec() : vec(std::vector<double>()) { }
        Phenovec(const std::vector<double> & v) : vec(v) { }
        virtual ~Phenovec() { }
        void push_back(const double el) { vec.push_back(el); }
        double operator[](unsigned int n) const { return vec[n]; }
        double& operator[](unsigned int n) { return vec[n]; }
        std::vector<double>::iterator begin() { return vec.begin(); }
        std::vector<double>::const_iterator begin() const { return vec.begin(); }
        std::vector<double>::iterator end() { return vec.end(); } 
        std::vector<double>::const_iterator end() const {return vec.end(); }
        unsigned int size() const { return vec.size(); }
        unsigned int dimensionality() const { return size(); }
    
    protected:
        std::vector<double> vec;
        
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & vec;
        }      
};


class Phenotype
{
	friend std::ostream& operator << (std::ostream&, const Phenotype&);
	friend void outformat(std::ostream &, const Phenotype &, 
			unsigned int width, unsigned int precision, std::string sep);
	friend class Environment;
	
	public:
		//constructors/destructor
		Phenotype();
		Phenotype(const double);
		Phenotype(const std::vector<double> &);
		Phenotype(const Phenovec &);
		Phenotype(const std::vector<double> &, const std::vector<double> &);
		Phenotype(const Phenotype &);
		~Phenotype();
		
		//initialization
		void initialize(const std::vector<double>&);
		
		//operator overload
		Phenotype& operator= (const Phenotype &);
		
		// getters
		double operator[] (const unsigned int index) const;
		Phenovec get_pheno() const;
		double get_pheno(const unsigned int index) const;
		Phenovec get_unstab() const;
		double get_unstab(const unsigned int index) const;
						
		// number of phenotypes
		unsigned int dimensionality() const;	    
		
	protected:
		Phenovec pheno;
		Phenovec unstabpheno;
        
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & pheno;
            ar & unstabpheno;
        }  
};


///////////// Specific class to get multivariate statistics on phenotypes. 
class PhenotypeStat
{
	public:
		//constructor
		PhenotypeStat() = delete;
		PhenotypeStat(const PhenotypeStat &) = delete;
		PhenotypeStat(const std::vector<Phenotype> &);
		PhenotypeStat(const std::vector<Phenovec> &);
		~PhenotypeStat();
		
		//functions
		Phenovec means_phen() const; 
		Phenovec vars_phen() const; 
		
		Phenovec means_unstab() const;
		Phenovec means_log_unstab() const;		
		Phenovec vars_unstab() const;
		Phenovec vars_log_unstab() const;
		
		unsigned int dimensionality() const;
		
		static std::vector<std::vector<double> > transpose_phen_matrix(const std::vector<Phenotype> &);
		static std::vector<std::vector<double> > transpose_unstabphen_matrix(const std::vector<Phenotype> &);
		static std::vector<std::vector<double> > transpose_phenovec_matrix(const std::vector<Phenovec> &);
		
	protected:
		MultivariateStat * pheno;
		MultivariateStat * unstab;
};

class BivariatePhenovecStat 
{
	// At one point it might be useful to make it multivariate (TODO)
	public:
		// constructors
		BivariatePhenovecStat() = delete;
		BivariatePhenovecStat(const std::vector<Phenovec>&, const std::vector<Phenovec>&); // This is only bivariate
		~BivariatePhenovecStat();
	
		// functions
		Phenovec means(unsigned int) const;
		Phenovec vars(unsigned int) const;
		Phenovec cov(unsigned int, unsigned int) const;
		
	protected:
		std::vector<MultivariateStat *> cache;
};

#endif // PHENOTYPE_H_INCLUDED
