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
		void copy(const Phenotype &);
		
		//functions
		double operator[] (const unsigned int index) const;
		unsigned int dimensionality() const;
		void write_debug (std::ostream&) const;	
		void write_simple (std::ostream&) const;	

		// Output
	    friend std::ostream& operator << (std::ostream&, const Phenotype&); // probably doesn't need to be friend

		
	protected:
		std::vector<double> pheno;
};

///////////// Specific class to get multivariate statistics on phenotypes. 
class PhenotypeStat: public MultivariateStat
{
	public:
		PhenotypeStat(const std::vector<Phenotype> &);
		
		Phenotype means_phen() const {return(Phenotype(means()));}
		Phenotype vars_phen() const {return(Phenotype(vars()));}
	protected:
		static std::vector<std::vector<double> > transpose_phen_matrix(const std::vector<Phenotype> &);
};


#endif
