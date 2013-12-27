#ifndef ENVIRONMENT_H_INCLUDED
#define ENVIRONMENT_H_INCLUDED


#include "Parameters.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>


class Environment
{
public:
    // constructors / destructor
    Environment(const ParameterSet&);

    // initialization / instance
    static Environment* instance;
    static void initialize(const ParameterSet&);

    // functions
    static double rand_effect();
    static double get_sd();

private:
    double sd;

};



#endif // ENVIRONMENT_H_INCLUDED
