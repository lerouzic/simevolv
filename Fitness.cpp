#include "Fitness.h"
#include "main.h"

using namespace std;


// constructors and destructor

Fitness::Fitness(const ParameterSet& param)
    : param(param)
    , type(param.getpar(FITNESS_TYPE)->GetInt())
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

    switch(instance->param.getpar(FITNESS_FLUCT)->GetInt())
    {
    case 0: // No fluctuations
        instance->strength = s1;
        instance->optimum = o1;
        break;
    case 1: // Smooth fluctuations
        instance->strength = s2+(s1-s2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;
        instance->optimum = o2+(o1-o2)*(1.0+cos(2.0*generation_number*M_PI/double(T)))/2.0;
        break;
    case 2: // Periodic flips
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
        break;
    case 3: // Stochastic flips
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
        break;
    case 4: // Brownian motion
        if ((generation_number % T) == 0)
        {
            instance->strength = instance->strength + Random::randgauss()*s2;
            instance->optimum = instance->optimum +  Random::randgauss()*o2;
        }
        break;
    case 5: // White noise
        if ((generation_number % T) == 0)
        {
            instance->strength = s1 + Random::randgauss()*s2;
            instance->optimum = o1 + Random::randgauss()*o2;
        }
        break;
    default:
        assert("Fluctuating selection type unknown.");
        break;
    }
}


void Fitness::update_extra(double strength)
{
    instance->strength = strength;
    instance->type = 1; // linear directional selection
}


double Fitness::compute(double phenotype, const Population & population)
{
    double population_value = Fitness::GetPopulationValue(population);
    return(compute(phenotype, population_value));
}


double Fitness::compute(double phenotype, double population_value)
{
    assert (instance != NULL);
    double fit = 0.0;

    //switch(instance->param.getpar(FITNESS_TYPE)->GetInt())
    switch(instance->type)
    {
    case 0 : // NoSel
        fit = 1.0;
        break;
    case 1 : // Linear
        fit = 1.0 + instance->strength*(phenotype - population_value);
        if (fit < 0)
        {
            fit = 0;
        }
        break;
    case 2 : // Exponential
        fit = exp(instance->strength*(phenotype - population_value));
        break;
    case 7 : // concave // take care, unordered!
        fit = 1.0 + 0.5*log(1.0+2.0*instance->strength*(phenotype - population_value));
        if (fit < 0.0)
            fit = 0.0;
        break;
    case 3 : // Gaussian
        fit = exp(-instance->strength*
                  (phenotype - instance->optimum)*
                  (phenotype - instance->optimum));
        break;
    case 4 : // Quadratic
        fit = 1.0 - instance->strength*
              (phenotype - instance->optimum)*
              (phenotype - instance->optimum);
        if (fit < 0)
        {
            fit = 0;
        }
        break;
    case 8 : // convex bilateral -- take care, unordered!
        fit = exp(-sqrt(instance->strength*
                        (phenotype - instance->optimum)*
                        (phenotype - instance->optimum)));
        break;
    case 5 : // Truncation Up
        fit = 0.0;
        if (phenotype > population_value)
            fit = 1.0;
        break;
    case 6 : // Truncation Down
        fit = 0.0;
        if (phenotype < population_value)
            fit = 1.0;
        break;
    default :
        assert("Fitness type unknown.");
        break;
    }
    return(fit);
}


double Fitness::GetPopulationValue(const Population& popul)
{
    vector<double> pheno = popul.phenotypes();

    switch(instance->param.getpar(FITNESS_TYPE)->GetInt())
    {
    case 1 :
    case 2 :
        return(popul.mean_phenotype());
        break;
    case 5 :
        sort(pheno.begin(), pheno.end());
        return(pheno[int(instance->strength*pheno.size())]);
        break;
    case 6 :
        sort(pheno.begin(), pheno.end());
        return(pheno[int((1-instance->strength)*pheno.size())]);
        break;
    default:
        break;
    }
    return(0);
}
