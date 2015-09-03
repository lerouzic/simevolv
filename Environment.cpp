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
    : sd_init(param.getpar(ENVIRO_SDINIT)->GetDouble())
    , sd_dynam(param.getpar(ENVIRO_SDDYNAM)->GetDouble())
    , sd_final(param.getpar(ENVIRO_SDFINAL)->GetDouble())
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

double Environment::init_disturb() 
{
	if (instance->sd_init == 0.)
		return(0.0);
	else
		return(instance->sd_init*Random::randgauss());
	
}

double Environment::dynam_disturb() 
{
	if (instance->sd_dynam == 0.)
		return(0.0);
	else
		return(instance->sd_dynam*Random::randgauss());
}

double Environment::final_disturb() 
{
	if (instance->sd_final == 0.)
		return(0.0);
	else
		return(instance->sd_final*Random::randgauss());
}
