#ifndef PHENOTYPE_H_INCLUDED
#define PHENOTYPE_H_INCLUDED

#include <vector>

class Phenotype
{
	public:
		Phenotype();
		Phenotype(const double);
		Phenotype(const std::vector<double> &);
		Phenotype(const Phenotype &);
		~Phenotype();
		
		Phenotype& operator= (const Phenotype &);
				
		double operator[] (const unsigned int index) const;
		unsigned int dimensionality() const;
		
		void write_debug (std::ostream&) const;	
		void write_simple (std::ostream&) const;	
		
	protected:
		void initialize(const std::vector<double>&);
		void copy(const Phenotype &);
	
		std::vector<double> pheno;
};



#endif
