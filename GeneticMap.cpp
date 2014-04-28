// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
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
#include "Architecture.h"
#include "Parconst.h"
#include "main.h"

#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

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


// functions

/* return the number of loci in the system */
int GeneticMap::nb_loc() const
{
    int nloc = (recrate.size()+1);

    return nloc;
}


/* calculate and return the recombination rate for the given loci */
double GeneticMap::recombination_rate(int loc1, int loc2) const
{
    int nloc = GeneticMap::nb_loc();

    if (loc2 == -1)
    {
        loc2 = loc1 + 1;
    }

    assert(loc1 >= 0);
    assert(loc2 < nloc);

    if (loc2 == loc1 + 1)
    {
        return(recrate[loc1]);
    }
    else
    {
        assert(1); // does not work yet (is it useful anyway?)
        double dist=0;
        for (int i = loc1; i < loc2-1; i++)
        {
            dist += recrate[loc1];
        }
        return(exp(-2*dist));
    }
}

