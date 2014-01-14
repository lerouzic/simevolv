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



#include "Fitness.h"
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

using namespace std;



// constructors and destructor

Fitness::Fitness(const ParameterSet& param)
    : param(param)
    , type(param.getpar(FITNESS_TYPE)->GetString())
    , strength(param.getpar(FITNESS_STRENGTH)->GetDouble())
    , optimum(param.getpar(FITNESS_OPTIMUM)->GetDouble())
{
}


// instance and initialization

Fitness * Fitness::instance = NULL;


void Fitness::initialize(const ParameterSet& param)
{
    if (Fitness::instance != NULL)
    {
        delete Fitness::instance;
        Fitness::instance = NULL;
    }
    Fitness::instance = new Fitness(param);
}


// functions

void Fitness::update_generation(const long unsigned int generation_number)
{
    long int T = instance->param.getpar(FITNESS_PERIOD)->GetInt();
    double s1 = instance->param.getpar(FITNESS_STRENGTH)->GetDouble();
    double s2 = instance->param.getpar(FITNESS_STRENGTH2)->GetDouble();
    double o1 = instance->param.getpar(FITNESS_OPTIMUM)->GetDouble();
    double o2 = instance->param.getpar(FITNESS_OPTIMUM2)->GetDouble();

	string fluct_type = instance->param.getpar(FITNESS_FLUCT)->GetString();
	if (fluct_type == FF_nofluct) 
	{
        instance->strength = s1;
        instance->optimum = o1;		
	} 
	else if (fluct_type == FF_smooth) 
	{
        instance->strength = s2+(s1-s2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;
        instance->optimum = o2+(o1-o2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;		
	} 
	else if (fluct_type == FF_pflips) 
	{
        if (int(generation_number / double(T/2)) % 2 == 0)
        {
            instance->strength = s1;
            instance->optimum  = o1;
        }
        else
        {
            instance->strength = s2;
            instance->optimum  = o2;
        }		
	} 
	else if (fluct_type == FF_sflips) 
	{
        // It's simpler to ignore the current values and to switch independently:
        if (Random::randnum() < 1.0 / double(T) / 2.0)
        {
            instance->strength = s1;
            instance->optimum  = o1;
        }
        if (Random::randnum() < 1.0 / double(T) / 2.0)
        {
            instance->strength = s2;
            instance->optimum  = o2;
        }		
	} 
	else if (fluct_type == FF_brown) 
	{
        if ((generation_number % T) == 0)
        {
            instance->strength = instance->strength + Random::randgauss()*s2;
            instance->optimum = instance->optimum +  Random::randgauss()*o2;
        }		
	} 
	else if (fluct_type == FF_white) 
	{
        if ((generation_number % T) == 0)
        {
            instance->strength = s1 + Random::randgauss()*s2;
            instance->optimum = o1 + Random::randgauss()*o2;
        }		
	} 
	else 
	{
        assert("Fluctuating selection type unknown.");		
	}
}


void Fitness::update_extra(double strength)
{
    instance->strength = strength;
    instance->type = 1; // linear directional selection
}


double Fitness::compute(const Phenotype & phenotype, const Population & population)
{
    double population_value = Fitness::GetPopulationValue(population);
    return(compute(phenotype, population_value));
}


double Fitness::compute(const Phenotype & phenotype, double population_value)
{
	int focal_phen = 0; // This will have to change for multivariate selection
	
    assert (instance != NULL);
    double fit = 0.0;

	string type = instance->type;
	if (type == FT_nosel) {
		fit = 1.0;
	} else if (type == FT_linear) {
        fit = 1.0 + instance->strength*(phenotype[focal_phen] - population_value);
        if (fit < 0)
        {
            fit = 0;
        }		
	} else if (type == FT_expo) {
		fit = exp(instance->strength*(phenotype[focal_phen] - population_value));
	} else if (type == FT_gauss) {
		fit = exp(-instance->strength*
			(phenotype[focal_phen] - instance->optimum)*
            (phenotype[focal_phen] - instance->optimum));		
	} else if (type == FT_quad) {
        fit = 1.0 - instance->strength*
              (phenotype[focal_phen] - instance->optimum)*
              (phenotype[focal_phen] - instance->optimum);
        if (fit < 0)
        {
            fit = 0;
        }		
	} else if (type == FT_truncup) {
        fit = 0.0;
        if (phenotype[focal_phen] > population_value)
            fit = 1.0;		
	} else if (type == FT_truncdown) {
        fit = 0.0;
        if (phenotype[focal_phen] < population_value)
            fit = 1.0;		
	} else if (type == FT_concave) {
        fit = 1.0 + 0.5*log(1.0+2.0*instance->strength*(phenotype[focal_phen] - population_value));
        if (fit < 0.0)
            fit = 0.0;		
	} else if (type == FT_convex) {
        fit = exp(-sqrt(instance->strength*
            (phenotype[focal_phen] - instance->optimum)*
            (phenotype[focal_phen] - instance->optimum)));		
	} else {
		assert("Fitness type unknown.");
	}

    return(fit);
}


double Fitness::GetPopulationValue(const Population& popul)
{
    //~ vector<double> pheno = popul.phenotypes();

	string type = instance->type;
	if (type == FT_linear || type == FT_expo) 
	{
		return(popul.mean_phenotype());
	} 
	else if (type == FT_truncup) 
	{
		// Broken
		assert("This option is broken");
        //~ sort(pheno.begin(), pheno.end());
        //~ return(pheno[int(instance->strength*pheno.size())]);
     } 
     else if (type == FT_truncdown) 
     {
		 assert("This option is broken");
		//~ sort(pheno.begin(), pheno.end());
        //~ return(pheno[int((1-instance->strength)*pheno.size())]);
     }
     return(0.0);
}
