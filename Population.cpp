// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "main.h"
#include "Population.h"
#include "Fitness.h"
#include "OutputFormat.h"
#include "Parconst.h"
#include "Statistics.h"
#include "Architecture.h"
#include "Canalization.h"
#include "Heritability.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

using namespace std;



// constructors and destructor

/* default constructor */
Population::Population()
{
}


/* constructor using the population size */
Population::Population(long int size)
{
    for (long int i = 0; i <= size; i++)
    {
        Individual indiv;
        pop.push_back(indiv);
    }
}


/* copy constructor */
Population::Population(const Population & copy)
    : pop(copy.pop)
    , nb_canal_test(copy.nb_canal_test)
    , nb_herit_test(copy.nb_herit_test)
{
}


/* constructor using a vector of individuals */
Population::Population(const std::vector<Individual>& vecindiv)
    : pop(vecindiv)
{
}


/* constructor using the parameters from the Paramers files */
Population::Population(const ParameterSet& param)
{
    initialize(param);
}



// operator overload

Population & Population::operator=(const Population& copy)
{
    if (this == &copy)
        return(*this);

    pop = copy.pop;
    nb_canal_test = copy.nb_canal_test; 
    nb_herit_test = copy.nb_herit_test;
    return(*this);
}


// instance and initialization

/* initialization of the population system and parameters */
void Population::initialize(const ParameterSet& param)
{	
    int popsize = param.getpar(INIT_PSIZE)->GetInt();
    //pop.resize(popsize);
    for (long int i = 0; i < popsize; i++)
    {
        Individual indiv(param);
        pop.push_back(indiv);
        //cout << i << endl;
    }
    nb_canal_test = param.getpar(OUT_CANAL_TESTS)->GetInt();
    nb_herit_test = param.getpar(OUT_HERIT_TESTS)->GetInt();
    update();
}


// functions

double fun_sqrt(double x) // I don't remember why this stupid function was necessary???
{
    return(std::sqrt(x));
}


/* reproduction of the population to give a new one */
Population Population::reproduce(long int offspr_number) const
{
    Population offspring;
    offspring.nb_canal_test = nb_canal_test; // Strange design, but otherwise the information is lost
    offspring.nb_herit_test = nb_herit_test;
    vector<double> cumul_fit = cumul_fitness();

    if (offspr_number == 0)
    {
        offspr_number = size();
    }

    // offspring.pop.resize(offspr_number);

    for (long int i = 0; i < offspr_number; i++)
    {
        offspring.pop.push_back(Individual::mate(
                               this->pick_parent(cumul_fit),
                               this->pick_parent(cumul_fit)));
    }
    offspring.update();
    return(offspring);
}


/* ???? */
void Population::update(void)
{
    double popvalue = Fitness::GetPopulationValue(*this);
    for (vector<Individual>::iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->update_fitness(popvalue);
    }
}


//~ vector<double> Population::phenotypes() const
//~ {
    //~ vector<double> pheno;
    //~ for (vector<Individual>::const_iterator indiv = pop.begin();
            //~ indiv != pop.end(); indiv++)
    //~ {
        //~ pheno.push_back(indiv->get_phenotype());
    //~ }
    //~ return(pheno);
//~ }


/* calculate and return the phenotypic mean value */
double Population::mean_phenotype() const
{
	int focal_phen = 0; // Dirty, needs to be fixed at one point
	vector<double> phen(pop.size());
	for (unsigned int i = 0; i < pop.size(); i++) {
		phen.push_back(pop[i].get_phenotype()[focal_phen]);
	}
	UnivariateStat us(phen);
    return(us.mean());
}


/* return the size of the population */
long int Population::size() const
{
    return(pop.size());
}


/* ?????? */ 
vector<double> Population::cumul_fitness() const
{
    vector<double> cum_fit(this->size());
    double cumul = 0;

    vector<double>::iterator cc = cum_fit.begin();
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end();
            indiv++, cc++)
    {
        cumul += indiv->get_fitness();
        *cc = cumul;
    }
    for (vector<double>::iterator i = cum_fit.begin();
            i != cum_fit.end(); i++)
    {
        *i = *i/cumul;
    }
    return(cum_fit);
}


/* ???? */
const Individual & Population::pick_parent(const vector<double>& cumfit) const
{
    // return(iterator_search_fit_table(rnum, cumfit));
    int i = search_fit_table(Random::randnum(), cumfit);
    return(pop[i]);
}


/* ???? */
long int Population::search_fit_table(double rnum, const vector<double>& cumfit) const
{
	long int i = stl_search_fit_table(rnum, cumfit);
	// this assertion takes time, and it seems to work: no need to check by default
	// assert(i == sequential_search_fit_table(rnum, cumfit));
    return(i);
}


/* ???? */
long int Population::sequential_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    long int i = 0;
    while (cumfit[i] < rnum)
        i++;
    return(i);
}


/* ???? */
long int Population::stl_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    vector<double>::const_iterator solution = std::lower_bound(cumfit.begin(), cumfit.end(), rnum);
    assert(solution != cumfit.end());
    return(solution - cumfit.begin());
}


/* determines if there will be a mutation in the population */
void Population::draw_mutation()
{
    for (unsigned int i = 0; i < pop.size(); i++) {
        pop[i].draw_mutation();
    }
    // population.pop.draw_mutation(population.pop);
}


/* force to make a mutation in the population */
void Population::make_mutation()
{
    int ind = floor(Random::randnum()*pop.size());
    pop[ind].make_mutation();
}


// output

void Population::write(int generation) const
{
    if (!OutputFormat::isInitialized())
    {
        cerr << "Warning: No output!\n";
    }
    write_debug (OutputFormat::GetDebug());
    write_xml   (OutputFormat::GetXml());
    write_simple(OutputFormat::GetSimple());
    write_summary(OutputFormat::GetSummary(), generation);
}


void Population::write_debug(ostream & out) const
{
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->write_debug(out);
    }
}


void Population::write_xml(ostream & out) const
{
    out << "xml output: not implemented yet.\n";
}


void Population::write_simple(ostream & out) const
{
    for (vector<Individual>::const_iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->write_simple(out);
    }
}


/* final summary of an analysis */
void Population::write_summary(ostream & out, int generation) const
{
	vector<Phenotype> phen;
	vector<Phenotype> gen;
	vector<double> fit;
	
	for (unsigned int i = 0; i < pop.size(); i++) {
		phen.push_back(pop[i].get_phenotype());
		gen.push_back(pop[i].get_genot_value());
		fit.push_back(pop[i].get_fitness());
	}
	
	PhenotypeStat phenstat(phen);
	PhenotypeStat genstat(gen);
	UnivariateStat fitstat(fit);	
	
	if (generation==1)
	{ 
		out << "Gen" << "\t";
		for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {
			out << "MPhen" << i + 1 << "\t";
		}
		for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {
			out << "VPhen" << i + 1 << "\t";
		}
		out << "MFit" << "\t";
		out << "VFit" << "\t";
		out << "FitOpt" << "\t";
		if (nb_canal_test > 0) 
		{
			for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {
				out << "CanPhen" << i + 1 << "\t";
			}
			out << "CanalFit" << "\t";
		}
		if (nb_herit_test > 0)
		{	
			for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {			
				out << "HerPhen" << i + 1 << "\t";
			}
			out << "HeritFit" << "\t";
		}
		out << endl; 
   	}
		
	out << generation << "\t";
    out << phenstat.means_phen() << "\t";
    out << phenstat.vars_phen() << "\t";
    out << fitstat.mean() << "\t";
    out << fitstat.var() << "\t";
    out << Fitness::current_optimum() << "\t";
    if (nb_canal_test > 0) {
		Canalization can_test(nb_canal_test, *this);
		out << can_test.phen_canalization() << "\t";
		out << can_test.fitness_canalization() << "\t";
	}    
	if (nb_herit_test > 0) {
		Heritability herit_test(nb_herit_test, *this);
		out << herit_test.h2() << "\t";
		out << herit_test.fit_h2() << "\t";
	}
    out << endl;
}

