// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

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


class Environment
{	
	// Singleton pattern
	
	public:
	    Environment(const ParameterSet&);
		~Environment() { }
		
	    // initialization / instance
	    // warning: dangerous structure, the instance is public! 
	    static Environment* instance;
	    static void initialize(const ParameterSet&);
	
	    // functions
		static double init_disturb(bool test = false);
		static double dynam_disturb(bool test = false);
		static double final_disturb();
	
	private:
	    double sd_init;
	    double sd_init_test;
	    double sd_dynam;
	    double sd_dynam_test;
	    double sd_final;	
};

#endif // ENVIRONMENT_H_INCLUDED
