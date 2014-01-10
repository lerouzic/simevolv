#ifndef CANALIZATION_H_INCLUDED
#define CANALIZATION_H_INCLUDED

#include "Phenotype.h"
#include "Statistics.h"
#include "Parameters.h"

#include <vector>

class Individual; // crossed references: Population.h -> Canalization.h -> Individual.h -> Population.h
class Population;

class Canalization
{
	public:
		Canalization();
		Canalization(const ParameterSet &); // Allows to store a meaningful value for nb_tests
		Canalization(unsigned int);
		~Canalization();
		
		unsigned int nb_tests() const {return _nb_tests;}
		
		// fill the object
		void reference_indiv(Individual);
		void mutant_indiv (Individual);
				
		// get the results
		Phenotype phen_canalization();
		double fitness_canalization();
				
	protected:
		void process();	
		void process_phen();
		void process_fit();
		
		std::vector<Individual> reference;
		std::vector<std::vector<Individual> > mutants;
		
		std::vector<Phenotype> mean_per_indiv;
		Phenotype mean_of_mean;
		Phenotype var_of_mean;
		
		std::vector<Phenotype> var_per_indiv;
		Phenotype mean_of_var;
		Phenotype var_of_var;
		
		std::vector<double> indiv_fitness_mean;
		std::vector<double> indiv_fitness_var;
		
		bool phen_ready;
		bool fit_ready;
		unsigned int _nb_tests;
	
};


#endif // CANALIZATION_H_INCLUDED
