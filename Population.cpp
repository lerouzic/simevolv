// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
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
#include "Individual.h"
#include "Fitness.h"
#include "OutputFormat.h"
#include "Parconst.h"
#include "Statistics.h"
#include "Architecture.h"
#include "Canalization.h"
#include "Heritability.h"
#include "Direpistasis.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>

using namespace std;



// constructors (note: most are useless)

/* default constructor (necessary a build a population individual by individual, but might reflect a design error) */ 
Population::Population()
	: nb_canal_test(0)
	, nb_herit_test(0)
	, nb_direpi_test(0)
{
}

/* copy constructor. Warning: it should be updated when the class changes, along with the assignement operator */ 
Population::Population(const Population & copy)
    : pop(copy.pop)
    , nb_canal_test(copy.nb_canal_test)
    , nb_herit_test(copy.nb_herit_test)
    , nb_direpi_test(copy.nb_direpi_test)
{
}


/* constructor using a vector of individuals */
Population::Population(const std::vector<Individual>& vecindiv)
    : pop(vecindiv)
    , nb_canal_test(0)
	, nb_herit_test(0)
	, nb_direpi_test(0)
{
}


/* constructor using the parameters from the Parameters files -- the only useful one, probably*/
Population::Population(const ParameterSet& param)
{
    initialize(param);
}



// assigment operator overload. Should be updated as the same time as the copy constructor. 

Population & Population::operator=(const Population& copy)
{
    if (this == &copy)
        return(*this);

    pop = copy.pop;
    nb_canal_test = copy.nb_canal_test; 
    nb_herit_test = copy.nb_herit_test;
    nb_direpi_test = copy.nb_direpi_test;
    return(*this);
}


// instance and initialization

/* initialization of the population from parameters */
void Population::initialize(const ParameterSet& param)
{	
    int popsize = param.getpar(INIT_PSIZE)->GetInt();
    for (long int i = 0; i < popsize; i++)
    {
        Individual indiv(param);
        pop.push_back(indiv);
    }
    nb_canal_test = param.getpar(OUT_CANAL_TESTS)->GetInt();
    nb_herit_test = param.getpar(OUT_HERIT_TESTS)->GetInt();
    nb_direpi_test = param.getpar(OUT_DIREPI_TESTS)->GetInt();
    update();
}


// functions

/* Sexual reproduction of the population 
     returns the offspring population (of size offspr_number) */
Population Population::reproduce(long int offspr_number /* = 0 */) const
{
    Population offspring;
    
    // Information about the number of canalization, heritability, and directionality tests
    // is transmitted to the offspring. Strange design, but the Population
    // class does not save a copy of the parameter set. 
    offspring.nb_canal_test = nb_canal_test; 
    offspring.nb_herit_test = nb_herit_test;
    offspring.nb_direpi_test = nb_direpi_test;
    
    // cumulated fitnesses. Computing it here fastens the random sampling algorithm.
    vector<double> cumul_fit = cumul_fitness();

    if (offspr_number == 0)
    {	// When the population size is not provided, it is expected to be
		// the same as in the parental population.
        offspr_number = this->size();
    }

    for (long int i = 0; i < offspr_number; i++)
    {
		// Each offspring results from a cross between two random parents.
		// (Hermaphrodite, sexual individuals)
        offspring.pop.push_back(Individual::mate(
                               this->pick_parent(cumul_fit),
                               this->pick_parent(cumul_fit)));
    }
    offspring.update();
    return(offspring);
}


/* Updates the status of the population. 
   So far, the only variables that change are individual fitnesses.
   For some fitness functions, the fitness can only be computed once the
   population is fully known */
/* Note: there would be a way to make this more elegant: disallow the
   default constructor, and initialize the population in one step from
   a vector of individuals. This may generate several minor issues, and
   is probably not necessary */
void Population::update(void)
{
    Fitness::update(*this);
	// Initialization of the fitness function for this generation
    for (vector<Individual>::iterator indiv = pop.begin();
            indiv != pop.end(); indiv++)
    {
        indiv->update_fitness(*this);
    }
}



/* calculates and returns the phenotypic mean value for each trait */
Phenovec Population::mean_phenotype() const
{
	vector<Phenotype> phen;
	for (unsigned int i = 0; i < pop.size(); i++) {
		phen.push_back(pop[i].get_phenotype());
	}
	PhenotypeStat phenstat(phen);
	return(phenstat.means_phen());
}


/* return the size of the population */
long int Population::size() const
{
    return(pop.size());
}


/* Generates a vector containing the cumulated sum of individual
   fitnesses, in the same order as in the pop vector.
   The purpose of this function is to fasten the reproduction algotrithm.
   The function scales fitnesses such as the sum is 1.0 */ 
vector<double> Population::cumul_fitness() const
{
    vector<double> cum_fit(this->size());
    double cumul = 0.0;

	// Warning, double loop over both individuals and cum_fit
	// This is a bit overcomplicated. Why using iterators at all? 
    auto cc = cum_fit.begin();
    for (auto indiv = pop.begin();
            indiv != pop.end();
            indiv++, cc++)
    {
        cumul += indiv->get_fitness();
        *cc = cumul;
    }
    for (auto i = cum_fit.begin();
            i != cum_fit.end(); i++)
    {
        *i = *i/cumul;
    }
    return(cum_fit);
}


/* Picks a parent randomly, proportionally to individual fitnesses. Requires a vector
   of cumulated fitnesses. This function is just a wrapper for search_fit_table
   (just in case several algorithms should be compared) */
const Individual & Population::pick_parent(const vector<double>& cumfit) const
{
    int i = search_fit_table(Random::randnum(), cumfit);
    return(pop[i]);
}


/* Returns the index of the individuals matching the random number 0 < rnum < 1, from the cumulated 
   fitness vector. Several algorithms can be used here. */
long int Population::search_fit_table(double rnum, const vector<double>& cumfit) const
{
	long int i = stl_search_fit_table(rnum, cumfit);
	// this assertion takes time, and it seems to work: no need to check by default
	// assert(i == sequential_search_fit_table(rnum, cumfit));
    return(i);
}


/* Old (and not efficient) search algorithm (sequential search: tries all sorted values until 
   finding the proper one.  */
long int Population::sequential_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    long int i = 0;
    while (cumfit[i] < rnum)
        i++;
    return(i);
}


/* New (and efficient) search algorithm based on a binary search (using the STL library). */
long int Population::stl_search_fit_table(double rnum, const vector<double>& cumfit) const
{
    auto solution = std::lower_bound(cumfit.begin(), cumfit.end(), rnum);
    assert(solution != cumfit.end());
    return(solution - cumfit.begin());  // Warning: this works with vector, but probably not with other STL containers. 
}


/* determines if there will be a mutation in the population */
/* (just calls draw_mutation() on every single individual, not very original) */
void Population::draw_mutation()
{
    for (unsigned int i = 0; i < pop.size(); i++) {
        pop[i].draw_mutation();
    }
}


/* forces a mutation in the population */
/* contrary to draw_mutation(), here there is no mutation rate involved: 
   a mutation MUST occur in a random individual */
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


/* Final summary of an analysis */
/* This is the function that will generate the output file */
/* Note that it is not only a display function: some potentially heavy calculation is run here */
void Population::write_summary(ostream & out, int generation) const
{
	vector<Phenotype> phen;
	vector<Phenotype> gen;
	vector<double> fit;
	
	// stores vectors of phenotypic values, genotypic values, and fitnesses. 
	for (unsigned int i = 0; i < pop.size(); i++) {
		phen.push_back(pop[i].get_phenotype());
		gen.push_back(pop[i].get_genot_value());
		fit.push_back(pop[i].get_fitness());
	}
	
	// Calls the multidimensional statistical routines for phenotypic and genotypic values
	// Calls the unidimensional routine for fitnesses. 
	PhenotypeStat phenstat(phen);
	PhenotypeStat genstat(gen);
	UnivariateStat fitstat(fit);	
	
	/* First generation: need to write the headers. 
	   This is not extremely clean. 
       * No warranty that the function is called all the time at the first
         generation. If not, there will be no headers.
       * This cannot be disabled, which can be annoying. 
       * There is no double check that the headers actually match with the
         content of the columns. This part of the code needs to be carefully
         synchronized with the rest of the function! */
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
		for (unsigned int opt = 0; opt < Fitness::current_optimum().size(); opt++) {
			out << "FitOpt" << opt << "\t";
		}
		if (nb_canal_test > 0) 
		{
			for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {
				out << "CanPhen" << i + 1 << "\t";
			}
			out << "CanFit" << "\t";
		}
		if (nb_herit_test > 0)
		{	
			for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {			
				out << "HerPhen" << i + 1 << "\t";
			}
			out << "HerFit" << "\t";
		}
		if (nb_direpi_test > 0) 
		{
			for (unsigned int i = 0; i < phenstat.dimensionality(); i++) {
				out << "DirPhen" << i + 1 << "\t";
			}
			out << "DirFit" << "\t";
		}
		out << endl; 
   	}
		
  /* The real results are written now. For a n-dimensional phenotype:
	 * Col 1: generation number
	 * n following columns: phenotypic means
	 * n following columns: phenotypic variances
	 * following column: fitness mean
	 * following column: fitness variance
	 * following column: fitness optimum. If not relevant, something meaningless will be written anyway
	 * n following columns: phenotype canalization score (if enabled)
	 * following column: fitness canalization score (if enabled)
	 * n following columns: heritability for the n phenotypes (if enabled)
	 * following column: heritability for fitness (if enabled)
  */
		
	out << generation << "\t";
	Phenovec mm = phenstat.means_phen();
	Phenovec mm2 = phenstat.means_unstab();
	for (unsigned int i = 0; i < mm.size(); i++){
		//out << mm[i] << "(" << mm2[i] << ")" << "\t";  //debug
		out << mm[i] << "\t";
	}
	Phenovec vv = phenstat.vars_phen();
	for (unsigned int i = 0; i < vv.size(); i++){
		out << vv[i] << "\t";
	}
    out << fitstat.mean() << "\t";
    out << fitstat.var() << "\t";
    out << Fitness::current_optimum() << "\t";
    if (nb_canal_test > 0) {
		// Runs the canalization tests
		Canalization can_test(nb_canal_test, *this);
		out << can_test.phen_canalization() << "\t";
		out << can_test.fitness_canalization() << "\t";
	}    
	if (nb_herit_test > 0) {
		// Runs the heritability tests
		Heritability herit_test(nb_herit_test, *this);
		out << herit_test.h2() << "\t";
		out << herit_test.fit_h2() << "\t";
	}
	if (nb_direpi_test > 0) {
		// Runs the directional epistasis tests
		Direpistasis dir_test(nb_direpi_test, *this);
		out << dir_test.phen_direpistasis() << "\t";
		out << dir_test.fitness_direpistasis() << "\t";
	}
    out << endl;
}

