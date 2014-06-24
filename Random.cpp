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



#include "Random.h"

#include <ctime>
#include <unistd.h>
#include <cassert>
#include <iostream>

using namespace std;



// constructors/destructor

/* default constructor */
Random::Random()
{
    random_generator = gsl_rng_alloc(gsl_rng_mt19937);
    seed = time(0)*getpid();
    //cerr << "Random number generator initialized with time():" << seed << endl << endl;
    gsl_rng_set (random_generator, seed);
}


/* constructor with seed (for repetability) */
Random::Random(long int s) : seed(s)
{
    random_generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (random_generator, seed);
}


/* destructor */
Random::~Random()
{
    gsl_rng_free(random_generator);
    random_generator = NULL;
    delete Random::instance;
}



// initialization

/* put the existence of the random number system to non-existent */
Random* Random::instance = NULL;


/* verify if the random system number is initialized */
bool Random::is_initialized()
{
    return(!(Random::instance == NULL));
}


/* initialize the random number system */
void Random::initialize()
{
    if (Random::instance != NULL)
    {
        delete Random::instance;
        Random::instance = NULL;
    }
    Random::instance = new Random();
}


/* initialize the random number system with a given seed */
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

/* get the seed */
long int Random::get_seed()
{
    assert(Random::instance != NULL);
    return(Random::instance->seed);
}


/* return a random number 
 * from a numeric distribution */
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


/* return a random number 
 * from a gaussian distribution */
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

