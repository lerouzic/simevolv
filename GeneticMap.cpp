// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "GeneticMap.h"
#include "Parconst.h"

#include <cassert>

using namespace std;


// constructors and destructor

/* default constructor */
GeneticMap::GeneticMap()
    : recrate(vector<double> (0))
{
}

/* constructor using the parameters from the Parameters files */
GeneticMap::GeneticMap(const ParameterSet& param)
    : recrate(vector<double> (0))
{
    int nloc = param.getpar(GENET_NBLOC) -> GetInt();

    for (int i = 0; i < nloc-1; i++)
    {
        recrate.push_back(param.getpar(GENET_RECRATES)->GetDouble(i));
    }
}

GeneticMap::~GeneticMap() 
{
}

// functions

/* return the number of loci in the system */
int GeneticMap::nb_loc() const
{
    return(recrate.size()+1);
}

/* calculate and return the recombination rate for the given locus */
double GeneticMap::recombination_rate(unsigned int loc) const
{
    unsigned int nloc = GeneticMap::nb_loc();
    assert(loc >= 0);
	assert(loc < nloc - 1);

	return(recrate[loc]);
}
