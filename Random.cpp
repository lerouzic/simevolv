#include "Random.h"

#include <ctime>
#include <unistd.h>
#include <cassert>
#include <iostream>

using namespace std;



// constructors/destructor

Random::Random()
{
    random_generator = gsl_rng_alloc(gsl_rng_mt19937);
    seed = time(0)*getpid();
    cerr << "Warning: Random number generator initialized with time():" << seed << endl << endl;
    gsl_rng_set (random_generator, seed);
}


Random::Random(long int s) : seed(s)
{
    random_generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (random_generator, seed);
}


Random::~Random()
{
    gsl_rng_free(random_generator);
    random_generator = NULL;
}



// initialization

Random* Random::instance = NULL;


bool Random::is_initialized()
{
    return(!(Random::instance == NULL));
}


void Random::initialize()
{
    if (Random::instance != NULL)
    {
        delete Random::instance;
        Random::instance = NULL;
    }
    Random::instance = new Random();
}


void Random::initialize(long int s)
{
    if (Random::instance != NULL)
    {
        delete Random::instance;
        Random::instance = NULL;
    }
    Random::instance = new Random(s);
}



// functions

long int Random::get_seed()
{
    assert(Random::instance != NULL);
    return(Random::instance->seed);
}


double Random::randnum()
{
    if (!Random::is_initialized())
    {
        cerr << "The random number generator is used prior to initialization" << endl;
    }
    if (Random::instance == NULL)
    {
        Random::initialize();
    }
    return(gsl_rng_uniform(Random::instance->random_generator));
}


double Random::randgauss()
{
    if (!Random::is_initialized())
    {
        cerr << "The random number generator is used prior to initialization" << endl;
    }
    if (Random::instance == NULL)
    {
        Random::initialize();
    }
    return(gsl_ran_gaussian(Random::instance->random_generator, 1.0));
}

