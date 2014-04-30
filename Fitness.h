// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef FITNESS_H_INCLUDED
#define FITNESS_H_INCLUDED

#include "Parameters.h"
#include "Population.h"
#include "Phenotype.h"

#include <string>


class Fitness
{
	// This is a singleton class: only one instance, access only through static
	// functions (no public constructor)
    public:
        // instance/initialization
        static void initialize(const ParameterSet&);

        // fonctions
        static double GetPopulationValue(const Population &);
        static void update_generation(const long unsigned int);
        static double compute(const Phenotype& phenotype, const Population & population);
        static double compute(const Phenotype& phenotype, double population_value);
        static double current_optimum() {return(instance->optimum);}

    protected:
        // constructors/destructor
        Fitness(const ParameterSet&);    
    
        static Fitness * instance;

        const ParameterSet & param;
        std::string type;
        double strength;
        double optimum;
};


#endif // FITNESS_H_INCLUDED

