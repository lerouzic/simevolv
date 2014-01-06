#ifndef PHENOTYPE_H_INCLUDED
#define PHENOTYPE_H_INCLUDED

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
		
		
	protected:
		std::vector<double> pheno;
};



#endif
