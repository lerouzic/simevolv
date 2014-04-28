// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



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
