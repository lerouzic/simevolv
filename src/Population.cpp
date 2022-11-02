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
#include "Random.h"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

using namespace std;

// constructors (note: most are useless)

/* default constructor (necessary a build a population individual by individual, but might reflect a design error) */ 
Population::Population()
	: selfing_rate(0.0)	
    , clonal_rate(0.0)
	, nb_canal_test(0)
	, nb_herit_test(0)
	, nb_direpi_test(0)
	, out_geno("no")
	, out_unstab("no")
{
}

Population::Population(const Population & copy)
    : pop(copy.pop)
	, selfing_rate(copy.selfing_rate)    
    , clonal_rate(copy.clonal_rate)
    , nb_canal_test(copy.nb_canal_test)
    , nb_herit_test(copy.nb_herit_test)
    , nb_direpi_test(copy.nb_direpi_test)
    , out_geno(copy.out_geno)
	, out_unstab(copy.out_unstab)
{
}

Population::Population(const std::vector<Individual>& vecindiv)
    : pop(vecindiv)
	, selfing_rate(0.0)    
    , clonal_rate(0.0)
    , nb_canal_test(0)
	, nb_herit_test(0)
	, nb_direpi_test(0)
	, out_geno("no")
	, out_unstab("no")
{
}

/* constructor using the parameters from the Parameters files -- the only useful one, probably*/
Population::Population(const ParameterSet& param)
{
    initialize(param);
}

Population & Population::operator=(const Population& copy)
{
    if (this == &copy)
        return(*this);

    pop = copy.pop;
 	selfing_rate = copy.selfing_rate;   
    clonal_rate = copy.clonal_rate;
    nb_canal_test = copy.nb_canal_test; 
    nb_herit_test = copy.nb_herit_test;
    nb_direpi_test = copy.nb_direpi_test;
    out_geno = copy.out_geno;
	out_unstab = copy.out_unstab;
    return(*this);
}


// instance and initialization

/* initialization of the population from parameters */
void Population::initialize(const ParameterSet& param)
{	
    int popsize = param.getpar(INIT_PSIZE)->GetInt();
    for (long int i = 0; i < popsize; i++)
    {
        Individual indiv(param); // Initialize the genotype only
        pop.push_back(indiv);
    }
    // Parameters that may change during simulation
	// update_param(param); // no harm, but not necessary here. 
    
    
	// Parameters that cannot change
    out_geno = param.getpar(OUT_GENO)->GetString(); 
    out_unstab = param.getpar(OUT_UNSTAB)->GetString();
}

/* Set up parameters. This may happen during initialization
 * or afterwards when the parameter file changes. */
void Population::update_param(const ParameterSet & param)
{
    selfing_rate = param.getpar(GENET_SELFING)->GetDouble();    
    clonal_rate  = param.getpar(GENET_CLONAL)->GetDouble();
    nb_canal_test = param.getpar(OUT_CANAL_TESTS)->GetInt();
    nb_herit_test = param.getpar(OUT_HERIT_TESTS)->GetInt();
    nb_direpi_test = param.getpar(OUT_DIREPI_TESTS)->GetInt();
    
    // Individuals don't know how to update themselves
    rate_type new_epigenet = param.getpar(GENET_EPIGENET)->GetDouble();
    for (auto &indiv : pop)
    {
        indiv.update_epigenet(new_epigenet);
    } 
}


void Population::update_phenotype(void)
{
    for (auto &indiv : pop)
    {
        indiv.update_phenotype();
    }    
}

void Population::update_fitness(void)
{
    Fitness::update(*this);
	// Initialization of the fitness function for this generation
    for (auto &indiv : pop)
    {
        indiv.update_fitness(*this);
    }
}

// functions

/* Sexual / asexual reproduction of the population 
     returns the offspring population (of size offspr_number) */
Population Population::reproduce(long int offspr_number /* = 0 */) const
{
    Population offspring;
    
    // Information about the number of canalization, heritability, and directionality tests
    // is transmitted to the offspring. Strange design, but the Population
    // class does not save a copy of the parameter set. 
    offspring.selfing_rate   = selfing_rate;
    offspring.clonal_rate    = clonal_rate; 
    offspring.nb_canal_test  = nb_canal_test; 
    offspring.nb_herit_test  = nb_herit_test;
    offspring.nb_direpi_test = nb_direpi_test;
    offspring.out_geno       = out_geno;
    offspring.out_unstab     = out_unstab;
    
    // cumulated fitnesses. Computing it here fastens the random sampling algorithm.
    vector<fitness_type> cumul_fit = cumul_fitness();

    if (offspr_number == 0)
    {	// When the population size is not provided, it is expected to be
		// the same as in the parental population.
        offspr_number = this->size();
    }

    for (long int i = 0; i < offspr_number; i++)
    {
		const Individual & firstpar = this->pick_parent(cumul_fit);
        		
        if (Random::randnum() < clonal_rate) { 
            // Clonal reproduction
            offspring.pop.push_back(firstpar.clone());
        } else { 
            // Sexual reproduction
            if (Random::randnum() < selfing_rate) {  
                // Self-fertilization: one parent
                offspring.pop.push_back(Individual::mate(firstpar, firstpar));
            } else { 
                // Outcrossing: two parents
                offspring.pop.push_back(Individual::mate(firstpar, this->pick_parent(cumul_fit)));
            }
        }
    }
    return(offspring);
}



/* calculates and returns the phenotypic mean value for each trait */
Phenotype Population::mean_phenotype() const
{
	vector<Phenotype> vphen;
	for (const auto & i : pop) 
		vphen.push_back(i.get_phenotype());

	return Phenotype::mean(vphen);
}

long int Population::size() const
{
    return(pop.size());
}

/* Generates a vector containing the cumulated sum of individual
   fitnesses, in the same order as in the pop vector.
   The purpose of this function is to fasten the reproduction algorithm.
   The function scales fitnesses such as the sum is 1.0 */ 
vector<fitness_type> Population::cumul_fitness() const
{
    vector<fitness_type> cum_fit(pop.size());
    fitness_type cumul = 0.0;

	for (unsigned int i = 0; i < pop.size(); i++) 
	{
		cumul += pop[i].get_fitness();
		cum_fit[i] = cumul;
	}
	
	for (unsigned int i = 0; i < pop.size(); i++)
	{
		cum_fit[i] /= cumul;
	}
    return(cum_fit);
}

/* Picks a parent randomly, proportionally to individual fitnesses. Requires a vector
   of cumulated fitnesses. This function is just a wrapper for search_fit_table
   (just in case several algorithms should be compared) */
const Individual & Population::pick_parent(const vector<fitness_type>& cumfit) const
{
    int i = search_fit_table(Random::randnum(), cumfit);
    return(pop[i]);
}

/* Returns the index of the individuals matching the random number 0 < rnum < 1, from the cumulated 
   fitness vector. Several algorithms can be used here. */
long int Population::search_fit_table(fitness_type rnum, const vector<fitness_type>& cumfit) const
{
    return(stl_search_fit_table(rnum, cumfit));
}

/* Old (and not efficient) search algorithm (sequential search: tries all sorted values until 
   finding the proper one.  */
long int Population::sequential_search_fit_table(fitness_type rnum, const vector<fitness_type>& cumfit) const
{
    long int i = 0;
    while (cumfit[i] < rnum)
        i++;
    return(i);
}

/* New (and efficient) search algorithm based on a binary search (using the STL library). */
long int Population::stl_search_fit_table(fitness_type rnum, const vector<fitness_type>& cumfit) const
{
    auto solution = std::lower_bound(cumfit.begin(), cumfit.end(), rnum);
    assert(solution != cumfit.end());
    return(solution - cumfit.begin());  // Warning: this works with vector, but probably not with other STL containers. 
}

/* determines if there will be a mutation in the population */
/* (just calls draw_mutation() on every single individual, not very original) */
/* Note: probably useless, as real mutations in simulations are drawn during
 * reproduction */
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

/* Final summary of an analysis */
/* This is the function that will generate the output file for the simulation */
/* Note that it is not only a display function: some potentially heavy calculation is run here */
void Population::write_summary(ostream & out, int generation) const
{
	vector<Phenotype> phen;
	vector<fitness_type> fit;
	vector<rate_type> epi;
	vector<Phenotype> genot; // Strange...
	    
	// stores vectors of phenotypic values, genotypic values, fitnesses. 
	for (const auto & i : pop) 
	{
		phen.push_back(i.get_phenotype());
		fit.push_back(i.get_fitness());
		epi.push_back(i.get_epigenet());
                
        // This is not clean. Stores genotype as a phenotype
		vector<pheno_type> matrix_vector_indiv;
        for (unsigned int loc = 0 ; loc < Architecture::Get()->nb_loc() ; loc++) 
        {
            vector<pheno_type> tmp = i.genotype->combine_at_loc(loc, &Allele::output_mean);

            for (auto j : tmp) // fills the matrix
                matrix_vector_indiv.push_back(j);

        }
       genot.push_back(Phenotype(matrix_vector_indiv)); 	
	}
			 
	// Calls the unidimensional routine for fitnesses. 
	UnivariateStat<fitness_type> fitstat(fit);
	UnivariateStat<rate_type> epistat(epi);
			
	/* First generation: need to write the headers. 
	   This is not extremely clean. 
       No warranty that the function is called all the time at the first generation. If not, there will be no headers.
       This cannot be disabled, which can be annoying. 
       There is no double check that the headers actually match with the content of the columns. This part of the code needs to be carefully synchronized with the rest of the function! */
       
    const size_t dim_phen = phen[0].dimensionality();   // Clearly not elegant
    const size_t dim_gen  = genot[0].dimensionality();
       
	if (generation==0)
	{ 
		outformat(out, "Gen");
		for (unsigned int i = 0; i < dim_phen; i++) 
		{
			outformat(out, i+1, "MPhen");
		}
		
		//Variance already calculated with the covariance
		
			for (unsigned int i = 0; i < dim_phen; i++) 
		{
			outformat(out, i+1, "VPhen");
		}
		
		for (unsigned int i = 0; i < dim_phen; i++)
		{
			for (unsigned int j = 0; j < dim_phen; j++)
			{
				outformat(out, i+1, j+1, "CovPhen");
			}
		}
		
	
		
		
		
		if (out_unstab == OU_yes) 
		{
			for (unsigned int i = 0; i < dim_phen; i++) 
			{
				outformat(out, i+1, "MUnstab");
			}
			for (unsigned int i = 0; i < dim_phen; i++) 
			{
				outformat(out, i+1, "VUnstab");
			}
		}
		outformat(out, "MFit");
		outformat(out, "VFit");
		for (unsigned int opt = 0; opt < Fitness::current_optimum().size(); opt++) 
		{
			outformat(out, opt+1, "FitOpt");
		}
		outformat(out, "MEpi");
		outformat(out, "VEpi");		
		if (nb_canal_test > 0) 
		{
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "MGenCan");
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "VGenCan");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "MGenCov");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "VGenCov");
			outformat(out, "GenCanFit");
			
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "MIniCan");
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "VIniCan");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "MIniCov");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "VIniCov");
			outformat(out, "IniCanFit");
						
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "MDynCan");
			for (unsigned int i = 0; i < dim_phen; i++) 
				outformat(out, i+1, "VDynCan");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "MDynCov");
			for (unsigned int i = 0; i < dim_phen; i++)
				for (unsigned int j = 0; j < dim_phen; j++)
					outformat(out, i+1, j+1, "VDynCov");
			outformat(out, "DynCanFit");
		}
		if (nb_herit_test > 0)
		{	
			for (unsigned int i = 0; i < dim_phen; i++) 
			{			
				outformat(out, i+1, "HerPhen");
			}
			outformat(out, "HerFit");
			for (unsigned int i = 0; i < dim_phen; i++) 
			{
				for (unsigned int j = 0; j < dim_phen; j++)
				{
					outformat(out, i+1, j+1, "POCov");
				}
			}
		}
		if (nb_direpi_test > 0) 
		{
			for (unsigned int i = 0; i < dim_phen; i++) 
			{
				outformat(out, i+1, "DirPhen");
			}
			outformat(out, "DirFit");
		}
		if (out_geno == OG_yes) 
		{
			for (unsigned int i = 0; i < dim_gen; i++)
			{
				outformat(out, i+1, "MeanAll");
			}
			for (unsigned int i = 0; i < dim_gen; i++)
			{
				outformat(out, i+1, "VarAll");
			}
		}
		
		
		
		out << endl; 
   	}
		
  /* The real results are written now. For a n-dimensional phenotype:
	 * Col 1: generation number
	 * n following columns: phenotypic means
	 * n following columns: phenotypic variances
	 * n following columns: unstability means (if enabled)
	 * n following columns: unstability variances (if enabled)
	 * following column: fitness mean
	 * following column: fitness variance
	 * following column: fitness optimum. If not relevant, something meaningless will be written anyway
	 * following column: epigenetic mean
	 * following column: epigenetic variance
	 * n following columns: canalization means(if enabled)
	 * n following columns: canalization variances (if enabled)
	 * following column: fitness canalization score (if enabled)
	 * n following columns: initdisturb canalization means(if enabled)
	 * n following columns: initdisturb canalization variances (if enabled)
	 * following column: fitness initdisturb canalization score (if enabled)
	 * n following columns: enviro canalization means(if enabled)
	 * n following columns: enviro canalization variances (if enabled)
	 * following column: fitness enviro canalization score (if enabled) 
	 * n following columns: heritability for the n phenotypes (if enabled)
	 * following column: heritability for fitness (if enabled)
	 * n following colums: w-matrix means
	 * n following colums: w-matrix variances
  */
	
	outformat(out, generation);
	Phenotype mm = Phenotype::mean(phen);

    outformat(out, mm);

	Phenotype vv = Phenotype::var(phen);

    outformat(out, vv); // Variance already calculated with the covariance
    
    for (unsigned int i = 0; i < dim_phen; i++)
		outformat(out, Phenotype::vcov(phen, i));
    

	if ((out_unstab == OU_yes) || (out_unstab == OU_log))
    // Warning: log or natural scales are not decided here any longer
	{
        outformat2(out, mm);
        outformat2(out, vv);
	}

    outformat(out, fitstat.mean());
    outformat(out, fitstat.var());
    outformat(out, Fitness::current_optimum());
    outformat(out, epistat.mean());
    outformat(out, epistat.var());
    if (nb_canal_test > 0) 
    {
		// Runs the canalization tests
		GeneticCanalization can_test(nb_canal_test, *this);
		outformat(out, can_test.meanpop_canphen());
		outformat(out, can_test.varpop_canphen());
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, can_test.meanpop_vcov(i));
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, can_test.varpop_vcov(i));
		outformat(out, can_test.meanpop_canlogfit());
		
		DisturbCanalization disturb_test(nb_canal_test, *this);
		outformat(out, disturb_test.meanpop_canphen());
		outformat(out, disturb_test.varpop_canphen());
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, disturb_test.meanpop_vcov(i));
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, disturb_test.varpop_vcov(i));
		outformat(out, disturb_test.meanpop_canlogfit());	

		EnviroCanalization enviro_test(nb_canal_test, *this);
		outformat(out, enviro_test.meanpop_canphen());
		outformat(out, enviro_test.varpop_canphen());
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, enviro_test.meanpop_vcov(i));
		for(unsigned int i = 0; i < dim_phen; i++)
			outformat(out, enviro_test.varpop_vcov(i));
		outformat(out, enviro_test.meanpop_canlogfit());	
	} 
	if (nb_herit_test > 0) 
	{
		// Runs the heritability tests
		Heritability herit_test(nb_herit_test, *this);
		outformat(out, herit_test.h2());
		outformat(out, herit_test.fit_h2());
		for (unsigned int i = 0; i < dim_phen; i++)
			outformat(out, herit_test.covOffPar(i));
	}
	if (nb_direpi_test > 0) 
	{
		// Runs the directional epistasis tests
		Direpistasis dir_test(nb_direpi_test, *this);
		outformat(out, dir_test.phen_direpistasis());
		outformat(out, dir_test.fitness_direpistasis());
	}
	if (out_geno == OG_yes) 
	{
		Phenotype wm = Phenotype::mean(genot);
		Phenotype wv = Phenotype::var(genot);

		outformat(out, wm);
		outformat(out, wv);
	}
	out << '\n';

}

/* Detailed output of the population at a specific generation */
void Population::write_details(ostream & out) const
{
			
	/* Need to write the headers. 
	   This is not extremely clean. 
       This cannot be disabled, which can be annoying. 
       There is no double check that the headers actually match with the content of the columns. This part of the code needs to be carefully synchronized with the rest of the function! 
       Column names are suboptimal, as most need to match those of the write_summary function. */
       
    const size_t dim_phen = pop[0].get_phenotype().dimensionality();   // Clearly not elegant (the next is even worse)
    const size_t dim_gen  = Architecture::Get()->nb_loc() * pop[0].genotype->combine_at_loc(0, &Allele::combine_mean).size();
       
	for (unsigned int i = 0; i < dim_phen; i++) 
	{
		outformat(out, i+1, "MPhen");
	}
	
	if (out_unstab == OU_yes) 
	{
		for (unsigned int i = 0; i < dim_phen; i++) 
		{
			outformat(out, i+1, "MUnstab");
		}
	}
	
	outformat(out, "MFit");
	
	outformat(out, "MEpi");

	if (nb_canal_test > 0) 
	{
		for (unsigned int i = 0; i < dim_phen; i++) 
			outformat(out, i+1, "MGenCan");
			
		outformat(out, "GenCanFit");
		
		for (unsigned int i = 0; i < dim_phen; i++) 
			outformat(out, i+1, "MIniCan");
			
		outformat(out, "IniCanFit");
		
		for (unsigned int i = 0; i < dim_phen; i++) 
			outformat(out, i+1, "MDynCan");
			
		outformat(out, "DynCanFit");
	}
	
	if (nb_direpi_test > 0) 
	{
		for (unsigned int i = 0; i < dim_phen; i++) 
		{
			outformat(out, i+1, "DirPhen");
		}
		outformat(out, "DirFit");
	}
	if (out_geno == OG_yes) 
	{
		for (unsigned int i = 0; i < dim_gen; i++)
		{
			outformat(out, i+1, "MeanAll");
		}
	}

		out << endl; 

		
  /* The real results are written now. For a n-dimensional phenotype:
	 * Col 1: generation number
	 * n following columns: phenotypic values
	 * n following columns: unstability (if enabled)
	 * following column: fitness 
	 * following column: epigenetic value
	 * n following columns: canalization score(if enabled)
	 * following column: fitness canalization score (if enabled)
	 * n following columns: initdisturb canalization score(if enabled)
	 * following column: fitness initdisturb canalization score (if enabled)
	 * n following columns: enviro canalization score(if enabled)
	 * following column: fitness enviro canalization score (if enabled) 
	 * n following columns: directional epistasis for all traits (if enabled)
	 * following column: directional epistasis for fitness (if enabled)
	 * n following colums: w-matrix
  */
	
	for (const auto & i : pop) 
	{
		outformat(out, i.get_phenotype());
		
		if ((out_unstab == OU_yes) || (out_unstab == OU_log))
		// Warning: log or natural scales are not decided here any longer
		{
			outformat2(out, i.get_phenotype());
		}
		
		outformat(out, i.get_fitness());
		outformat(out, i.get_epigenet());
		if (nb_canal_test > 0) 
		{
			// Runs the canalization tests
			GeneticCanalization can_test(nb_canal_test, i, *this);
			outformat(out, can_test.meanpop_canphen());
			outformat(out, can_test.meanpop_canlogfit());

			DisturbCanalization disturb_test(nb_canal_test, i, *this);
			outformat(out, disturb_test.meanpop_canphen());
			outformat(out, disturb_test.meanpop_canlogfit());

			EnviroCanalization enviro_test(nb_canal_test, i, *this);
			outformat(out, enviro_test.meanpop_canphen());
			outformat(out, enviro_test.meanpop_canlogfit());
		} 

		if (nb_direpi_test > 0) 
		{
			// Runs the directional epistasis tests
			Direpistasis dir_test(nb_direpi_test, i, *this);
			outformat(out, dir_test.phen_direpistasis());
			outformat(out, dir_test.fitness_direpistasis());
		}
		
		if (out_geno == OG_yes) 
		{
			vector<pheno_type> matrix_vector_indiv;
			for (unsigned int loc = 0 ; loc < Architecture::Get()->nb_loc() ; loc++) 
			{
				vector<pheno_type> tmp = i.genotype->combine_at_loc(loc, &Allele::combine_mean);

				for (auto j : tmp) // fills the matrix
					matrix_vector_indiv.push_back(j);

			}
			outformat(out, Phenotype(matrix_vector_indiv)); 
		}
		out << '\n';
	}
	out << '\n';
}

#ifdef SERIALIZATION_TEXT
ostream & operator << (ostream & out, const Population & pop) {
        boost::archive::text_oarchive oa(out);
        oa << pop;
        return out;
}

istream & operator >> (istream & in, Population & pop) {
        boost::archive::text_iarchive ia(in);
        pop.pop.clear(); // not completely sure why this is necessary (is Individual default_constructible?)
        ia >> pop;
        return in;
}
#endif
