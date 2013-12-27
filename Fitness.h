#ifndef FITNESS_H_INCLUDED
#define FITNESS_H_INCLUDED

#include "Parameters.h"
#include "Random.h"
#include "Population.h"
#include "Phenotype.h"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>


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

