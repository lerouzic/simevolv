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

#include "types.h"
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
		static pheno_type init_disturb(bool test = false);
		static pheno_type dynam_disturb(bool test = false);
		static pheno_type final_disturb();
	
	private:
	    pheno_type sd_init;
	    pheno_type sd_init_test;
	    pheno_type sd_dynam;
	    pheno_type sd_dynam_test;
	    pheno_type sd_final;	
};

#endif // ENVIRONMENT_H_INCLUDED
