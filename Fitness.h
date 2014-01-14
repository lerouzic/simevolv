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



#ifndef FITNESS_H_INCLUDED
#define FITNESS_H_INCLUDED

#include "Parameters.h"
#include "Population.h"
#include "Phenotype.h"

#include <string>



class Fitness
{
    public:
        // constructors/destructor
        Fitness(const ParameterSet&);

        // instance/initialization
        static Fitness * instance;
        static void initialize(const ParameterSet&);

        // fonctions
        static void update_generation(const long unsigned int);
        static void update_extra(double strength);
        static double compute(const Phenotype& phenotype, const Population & population);
        static double compute(const Phenotype& phenotype, double population_value);
        static double GetPopulationValue(const Population &);
        static double current_optimum() {return(instance->optimum);}

    protected:
        const ParameterSet & param;
        std::string type;
        double strength;
        double optimum;
};


#endif // FITNESS_H_INCLUDED

