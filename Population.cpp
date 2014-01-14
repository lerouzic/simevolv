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



#include "main.h"
#include "Population.h"
#include "Fitness.h"
#include "OutputFormat.h"
#include "Parconst.h"
#include "Statistics.h"
#include "Architecture.h"
#include "Canalization.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

using namespace std;



// constructors and destructor

Population::Population()
{
}


Population::Population(long int size)
{
    for (long int i = 0; i <= size; i++)
    {
        Individual indiv;
        pop.push_back(indiv);
    }
}


Population::Population(const Population & copy)
    : pop(copy.pop)
    , canal_test(copy.canal_test.nb_tests())
{
}


Population::Population(const std::vector<Individual>& vecindiv)
    : pop(vecindiv)
{
}


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
    canal_test = copy.canal_test; // Not sure this is particularly clever...
    return(*this);
}


// instance and initialization

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
    canal_test = Canalization(param);
    update();
}


// functions

double fun_sqrt(double x) // I don't remember why this stupid function was necessary???
{
    return(std::sqrt(x));
}


Population Population::reproduce(long int offspr_number) const
{
    Population offspring;
    Canalization tmp_canal(canal_test.nb_tests());
    offspring.canal_test = tmp_canal; // Bad design, this kind of things should not happen
    vector<double> cumul_fit = cumul_fitness();

    if (offspr_number == 0)
    {
        offspr_number = size();
    }

    offspring.pop.resize(offspr_number);

    for (long int i = 0; i < offspr_number; i++)
    {
        offspring.pop[i] = Individual::mate(
                               this->pick_parent(cumul_fit),
                               this->pick_parent(cumul_fit));
    }
    offspring.update();
    return(offspring);
}


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


long int Population::size() const
{
    return(pop.size());
}


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


const Individual & Population::pick_parent(const vector<double>& cumfit) const
{
    // return(iterator_search_fit_table(rnum, cumfit));
    return(pop[search_fit_table(Random::randnum(), cumfit)]);
}


long int Population::search_fit_table(double rnum, const vector<double>& cumfit) const
{
    return(sequential_search_fit_table(rnum, cumfit));
}


long int Population::sequential_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    long int i = 0;
    while (cumfit[i] < rnum)
        i++;
    return(i);
}


Individual Population::iterator_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    // Does not work!
    // return(*(std::find_if(cumfit.begin(), cumfit.end(), std::bind2nd(std::less<double>(), rnum))));
    return(Individual());
}


void Population::draw_mutation()
{
    for (unsigned int i = 0; i < pop.size(); i++) {
        pop[i].draw_mutation();
    }
    // population.pop.draw_mutation(population.pop);
}


void Population::make_mutation()
{
    int ind = floor(Random::randnum()*pop.size());
    pop[ind].make_mutation();
}

void Population::canalization_test() const
// canal_test is mutable, it can thus be updated in this const function
{
	unsigned int nb_tests = canal_test.nb_tests();
	if (nb_tests > 0) {
		for (unsigned int i = 0; i < pop.size(); i++) {
			const Individual & ref = pop[i];
			canal_test.reference_indiv(ref);
			for (unsigned int test = 0; test < nb_tests; test++) {
					// the test needs to know the population (*this) to compute the fitness.
					// not very elegant, but I don't know how to do otherwise for selection
					// regimes that depend on the population
				canal_test.mutant_indiv(ref.test_canalization(1, *this));
			}
		} 
	} 
	// this function does not return anything, it just fills the object canal_test. 
}

// output

void Population::write() const
{
    if (!OutputFormat::isInitialized())
    {
        cerr << "Warning: No output!\n";
    }
    write_debug (OutputFormat::GetDebug());
    write_xml   (OutputFormat::GetXml());
    write_simple(OutputFormat::GetSimple());
    write_summary(OutputFormat::GetSummary());
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


void Population::write_summary(ostream & out) const
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
			
    out << "MeanPhen= " << phenstat.means_phen() << "\t";
    out << "VarPhen= " << phenstat.vars_phen() << "\t";
    out << "MeanFit= " << fitstat.mean() << "\t";
    out << "VarFit= " << fitstat.var() << "\t";
    out << "FitOpt= " << Fitness::current_optimum() << "\t";
    if (canal_test.nb_tests() > 0) {
		canalization_test();
		out << "Canalphen=" << canal_test.phen_canalization() << "\t";
		out << "Canalfit=" << canal_test.fitness_canalization() << "\t";
	}
    
    out << endl;
}

