#ifndef ARCHITECTURE_H_INCLUDED
#define ARCHITECTURE_H_INCLUDED

#include "Parameters.h"
#include "GeneticMap.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"

#include <iostream>
#include <vector>



class Architecture
{
	public :
	    //constructors/destructor
	    Architecture();
	    Architecture(const Architecture&);
	    Architecture (const ParameterSet&);
	    virtual ~Architecture(){};
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    // instance / initialization
	    static Architecture* instance;
	    static void initialize(const ParameterSet&);
	    static Architecture* Get();
	    static Architecture* Get(const ParameterSet*);
	    static Architecture* Get(const ParameterSet&);
	
	    //functions
	    int nb_loc() const;
	    int all_size() const;
	    double mutation_rate(int) const;
	    double mutation_sd(int) const;
	    double recombination_rate(int) const;
	    void draw_mutation(const Haplotype&) const;
	    double make_mutation(int, std::vector<Allele>) const;
		
		//inheritance
	    virtual double phenotypic_value(const Genotype&) const;
	
	protected :
	    GeneticMap gmap;
	    int nloc;
	    int sall;
	    std::vector<double> mutrate;
	    std::vector<double> mutsd;

};


#endif // ARCHITECTURE_H_INCLUDED

