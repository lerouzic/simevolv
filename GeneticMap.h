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



#ifndef GENETICMAP_H_INCLUDED
#define GENETICMAP_H_INCLUDED

#include "Parameters.h"

#include <iostream>
#include <vector>



class GeneticMap 
{
	public:
	    //constructors/destructor
	    GeneticMap();
	    GeneticMap(const ParameterSet&);
	
	    //functions
	    int nb_loc() const;
	    double recombination_rate(int loc1, int loc2 = -1) const;
	
	protected:
	    std::vector<double> recrate;
};


#endif // GENETICMAP_H_INCLUDED

