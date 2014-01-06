#include "Environment.h"
#include "Parconst.h"
#include "main.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

using namespace std;



// constructors and destructor

Environment::Environment(const ParameterSet & param)
    : sd(param.getpar(ENVIRO_SD)->GetDouble())
{
}


// instance and initialization

Environment * Environment::instance = NULL;


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

Phenotype Environment::rand_effect(const Phenotype& genot_values)
{ // Probably not very efficient, but this will probably not be used for complex architectures
  // in which the environmental effect is correlated to the genotype	 
    assert (Environment::instance != NULL);
    vector<double> updated;
    for (unsigned int i = 0; i < genot_values.dimensionality(); i++) {
		updated.push_back(genot_values[i]+Environment::instance->sd*Random::randgauss());
	}
	Phenotype newphen(updated);
    return(newphen);
}


double Environment::get_sd()
{
    assert (Environment::instance != NULL);
    return(Environment::instance->sd);
}
