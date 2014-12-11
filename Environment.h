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



#ifndef ENVIRONMENT_H_INCLUDED
#define ENVIRONMENT_H_INCLUDED

#include "Parameters.h"
#include "Phenotype.h"

// Singleton pattern
class Environment
{	
	public:
	    Environment(const ParameterSet&);
	
	    // initialization / instance
	    // warning: dangerous structure, the instance is public! 
	    static Environment* instance;
	    static void initialize(const ParameterSet&);
	
	    // functions
	    static Phenotype rand_effect(Phenotype);
	    static double get_sd();
	
	private:
	    double sd;
	
};


#endif // ENVIRONMENT_H_INCLUDED
