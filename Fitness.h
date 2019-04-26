// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2017 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef FITNESS_H_INCLUDED
#define FITNESS_H_INCLUDED

#include "types.h"
#include "Parameters.h"
#include "Phenotype.h"

#include <cassert>
#include <vector>
#include <gsl/gsl_matrix.h>


class Population;
// These classes are implemented afterwards for readability
class Fitness_Phenotype;
class Fitness_Stability;

typedef std::vector<pheno_type> FitnessOptimum;
typedef std::vector<fitness_type> FitnessStrength;

class Fitness
{
	// This is a singleton class: only one instance, access only through static
	// functions (no public constructor)
    public:
		Fitness(const Fitness &) = delete;
        static void Terminate();
        // instance/initialization
        static void initialize(const ParameterSet&);

        // fonctions
        static void update(const Population &);
        static fitness_type compute(const Phenotype&, const Population&);
        static FitnessOptimum current_optimum();
        		
    protected:
        // constructors/destructor
        Fitness(const ParameterSet&);
		~Fitness();
    
        static Fitness * instance;
        
        Fitness_Phenotype * fitphen;
        Fitness_Stability * fitstab;
};

/* Trivial tool to expand incomplete vectors. Here is the summary of the problem:
   For simplicity, a vector can be specific in the parameter file
   as a single number. In this case, all members of the vector would be considered as 
   identical. 
   Fitness strengths and optima are in such a situation. The problem is that there is 
   no easy way to know the size of the vector before the fitness functions are called in
   a real context. So, when the different fitness classes are initialized, vectors are stored
   as provided by the parameter set (either vectors or scalars). When used for the first time, 
   the size is compared to the dimensionality of phenotypes, and adjusted there -- this function
   is used for that purpose. 
   Obviously, this is not very elegant, but there is no reason that it shoudn't work. */

template<typename T>
std::vector<T> expand_vec(const std::vector<T>& templ, unsigned int maxsize) {
	assert(templ.size() > 0);
	std::vector<T> answer = templ;
	while(answer.size() < maxsize) 
	{
		answer.push_back(templ[0]);
	}
	return(answer);    
}

/*************************** Fitness_Phenotype ****************************/
/**                         Inheritance Scheme [(A): Abstract class]:
   Fitness_Phenotype (A)
      -> Fitness_Phenotype_Noselection
      -> Fitness_Phenotype_Directional (A)
          -> Fitness_Phenotype_Linear
          -> Fitness_Phenotype_Expo
          -> Fitness_Phenotype_Concave
      -> Fitness_Phenotype_Stabilizing (A)
          -> Fitness_Phenotype_Gaussian
          -> Fitness_Phenotype_Quadratic
          -> Fitness_Phenotype_Biconvex

**/

class Fitness_Phenotype
{ // Abstract class (interface)
	public:
		Fitness_Phenotype(const ParameterSet & param) { }
		virtual ~Fitness_Phenotype() = 0;
	
		virtual fitness_type get_fitness(const Phenotype&, const Population &);
		virtual fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&) = 0;
		virtual void update(const Population &) { }
		virtual FitnessOptimum get_optimum() const;
};

class Fitness_Phenotype_Noselection: public Fitness_Phenotype
{
	public:
		Fitness_Phenotype_Noselection(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&) { return(1.0); }
};

class Fitness_Phenotype_Directional: public Fitness_Phenotype
{ // Abstract class for directional selection
	public: 
		Fitness_Phenotype_Directional(const ParameterSet &);
		virtual ~Fitness_Phenotype_Directional() = 0;
		virtual fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population &) = 0;
		void update(const Population &);
	protected:
		FitnessStrength strength;
		Phenotype popmean;
		const Population * popmem;
};

class Fitness_Phenotype_Linear: public Fitness_Phenotype_Directional
{
	public: 
		Fitness_Phenotype_Linear(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Expo: public Fitness_Phenotype_Directional
{
	public: 
		Fitness_Phenotype_Expo(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Concave: public Fitness_Phenotype_Directional
{
	public:
		Fitness_Phenotype_Concave(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Stabilizing: public Fitness_Phenotype
{ // Abstract class for stabilizing selection
	public:
		Fitness_Phenotype_Stabilizing(const ParameterSet &);
		virtual ~Fitness_Phenotype_Stabilizing() { }
		virtual fitness_type get_fitness_trait(unsigned int, const Phenotype&, const Population&) = 0;
		FitnessOptimum get_optimum() const { return(optimum);}
				
	protected:
		FitnessStrength strength;
		FitnessOptimum optimum;
};

class Fitness_Phenotype_Gaussian: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Gaussian(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype&, const Population&);
};

class Fitness_Phenotype_MultivarGaussian: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_MultivarGaussian(const ParameterSet &);
		~Fitness_Phenotype_MultivarGaussian();
		fitness_type get_fitness(const Phenotype&, const Population &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype&, const Population&) { return (0.0); } // should not be called. 
		
	protected:
		void compute_invsigma(const size_t);
		
		FitnessStrength cor;      // vector of size (n^2-n)/2 containing correlation components. 
		gsl_matrix * invsigma;    // inverse of the variance-covariance matrix, should be stored. 
};

class Fitness_Phenotype_Quadratic: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Quadratic(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Biconvex: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Biconvex(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype&, const Population&);
};



/**************************** Fitness_Stability *********************************/
/**                           inheritance scheme:
      Fitness_Stability
          -> Fitness_Stability_Noselection
          -> Fitness_Stability_Exponential
**/

class Fitness_Stability
{ // Abstract class (interface)
	public: 
		Fitness_Stability(const ParameterSet &) { }
		virtual ~Fitness_Stability() { }
		
		fitness_type get_fitness(const Phenotype &);
		virtual fitness_type get_fitness_trait(unsigned int, const Phenotype&) = 0;
};

class Fitness_Stability_Noselection: public Fitness_Stability
{
	public:
		Fitness_Stability_Noselection(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype &);	
};

class Fitness_Stability_Exponential: public Fitness_Stability
{
	public:
		Fitness_Stability_Exponential(const ParameterSet &);
		fitness_type get_fitness_trait(unsigned int, const Phenotype &);
	protected:
		FitnessStrength strength;
};

#endif // FITNESS_H_INCLUDED

