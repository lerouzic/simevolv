#ifndef ENVIRONMENT_H_INCLUDED
#define ENVIRONMENT_H_INCLUDED

#include "Parameters.h"
#include "Phenotype.h"



class Environment
{
	public:
	    // constructors / destructor
	    Environment(const ParameterSet&);
	
	    // initialization / instance
	    static Environment* instance;
	    static void initialize(const ParameterSet&);
	
	    // functions
	    static Phenotype rand_effect(const Phenotype &);
	    static double get_sd();
	
	private:
	    double sd;
	
};



#endif // ENVIRONMENT_H_INCLUDED
