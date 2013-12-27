#include "Environment.h"
#include "Parconst.h"
#include "main.h"

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

double Environment::rand_effect()
{
    assert (Environment::instance != NULL);
    return(Environment::instance->sd*Random::randgauss());
}


double Environment::get_sd()
{
    assert (Environment::instance != NULL);
    return(Environment::instance->sd);
}
