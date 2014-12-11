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



#include "Environment.h"
#include "Parconst.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;



// constructors and destructor

/* constructor using the parameters from the parameters files */
Environment::Environment(const ParameterSet & param)
    : sd(param.getpar(ENVIRO_SD)->GetDouble())
{
}

// instance and initialization

/* put the existence of the environmental values to non-existent */
Environment * Environment::instance = NULL;


/* initialization of the environmental values system */
void Environment::initialize(const ParameterSet & param)
{
    if (instance != NULL)
    {
        delete Environment::instance;
        Environment::instance = NULL;
    }
    Environment::instance = new Environment(param);
}

// functions

/* modify the genotypic value by adding an environmental effect 
 * (Probably not very efficient, but this will probably not be used for complex 
 * architectures in which the environmental effect is correlated to the genotype) */
Phenotype Environment::rand_effect(Phenotype genot_values)
{	 
    assert (Environment::instance != NULL);
    for (unsigned int i = 0; i < genot_values.dimensionality(); i++) {
		genot_values.add_to_pheno(i, Environment::instance->sd*Random::randgauss());
	}
    return(genot_values);
}

/* return the effect of environmental modification */
double Environment::get_sd()
{
    assert (Environment::instance != NULL);
    return(Environment::instance->sd);
}
