#ifndef RANDOM_H
#define RANDOM_H

// GNU Scientific Library functions:
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



class Random
{
	public :
	    // constructors/destructor
	    Random();
	    Random(long int);
	    ~Random();
	
	    // initialization
	    static Random * instance;
	    static bool is_initialized();
	    static void initialize();
	    static void initialize(long int s);
	
	    // functions
	    static long int get_seed();
	    static double randnum();
	    static double randgauss();
	
	
	protected :
	    long int seed;
	    gsl_rng * random_generator;

};

#endif
