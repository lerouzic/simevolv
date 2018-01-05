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

#include "Parameters.h"
#include "Phenotype.h"

#include <cassert>
#include <vector>


class Population;
// These classes are implemented afterwards for readability
class Fitness_Fluct;
class Fitness_Phenotype;
class Fitness_Stability;

typedef double fitness_type;
typedef std::vector<pheno_type> FitnessOptimum;
typedef std::vector<fitness_type> FitnessStrength;

class Fitness
{
	// This is a singleton class: only one instance, access only through static
	// functions (no public constructor)
    public:
		Fitness(const Fitness &) = delete;
		~Fitness();
        // instance/initialization
        static void initialize(const ParameterSet&);

        // fonctions
        static void fluctuate(unsigned int);
        static void update(const Population &);
        static fitness_type compute(const Phenotype&, const Population&);
        static FitnessOptimum current_optimum();
        		
    protected:
        // constructors/destructor
        Fitness(const ParameterSet&);    
    
        static Fitness * instance;
        
        Fitness_Fluct * fitfluct;
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

/*************************** Fitness_Fluct ********************************/
/**                         Inheritance scheme [(A): abstract class]
    Fitness_Fluct (A)
        -> Fitness_Fluct_Nofluct
        -> Fitness_Fluct_States (A)
             -> Fitness_Fluct_Flips
             -> Fitness_Fluct_Sflips
             -> Fitness_Fluct_Smooth
        -> Fitness_Fluct_Noise (A)
             -> Fitness_Fluct_Whitenoise
             -> Fitness_Fluct_Brownian 
**/

class Fitness_Fluct
{ // Abstract class (interface)
	public:
		Fitness_Fluct(const ParameterSet & param) { }
		virtual ~Fitness_Fluct() = 0;
		
		virtual FitnessStrength get_new_strength(const FitnessStrength & old_strength, unsigned int generation) = 0;
		virtual FitnessOptimum get_new_optimum(const FitnessOptimum & old_optimum, unsigned int generation) = 0;
};

class Fitness_Fluct_Nofluct: public Fitness_Fluct
{
	public:
		Fitness_Fluct_Nofluct(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);		
};

class Fitness_Fluct_States: public Fitness_Fluct
{ // Abstract class
	public:
		Fitness_Fluct_States(const ParameterSet &);
		virtual ~Fitness_Fluct_States() = 0;  
		
	protected:
		FitnessStrength strength_state1;
		FitnessStrength strength_state2;
		FitnessOptimum optimum_state1;
		FitnessOptimum optimum_state2;
		unsigned int period;		
};

class Fitness_Fluct_Flips: public Fitness_Fluct_States
{
	public:
		Fitness_Fluct_Flips(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);
};

class Fitness_Fluct_Sflips: public Fitness_Fluct_States
{
	public: 
		Fitness_Fluct_Sflips(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);		
};

class Fitness_Fluct_Smooth: public Fitness_Fluct_States
{
	public:
		Fitness_Fluct_Smooth(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);	
};

class Fitness_Fluct_Noise: public Fitness_Fluct
{ // Abstract class
	public:
		Fitness_Fluct_Noise(const ParameterSet &);
		virtual ~Fitness_Fluct_Noise() = 0;	
			
	protected:
		FitnessStrength strength_sd;
		FitnessOptimum optimum_sd;
		unsigned int period;
};

class Fitness_Fluct_Whitenoise: public Fitness_Fluct_Noise
{
	public:
		Fitness_Fluct_Whitenoise(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);		
	protected:
		FitnessStrength strength_ref;
		FitnessOptimum optimum_ref;
};

class Fitness_Fluct_Brownian: public Fitness_Fluct_Noise
{
	public:
		Fitness_Fluct_Brownian(const ParameterSet &);
		FitnessStrength get_new_strength(const FitnessStrength &, unsigned int);
		FitnessOptimum get_new_optimum(const FitnessOptimum &, unsigned int);
};



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
	
		fitness_type get_fitness(const Phenotype&, const Population &);
		virtual fitness_type get_fitness_trait(unsigned int trait, const Phenotype&, const Population&) = 0;
		virtual void update(const Population &) { }
		virtual void fluctuate(Fitness_Fluct *, unsigned int) { }
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
		void fluctuate(Fitness_Fluct *, unsigned int);
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
		void fluctuate(Fitness_Fluct *, unsigned int);
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

