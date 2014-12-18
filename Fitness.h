// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

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
#include "Population.h"
#include "Phenotype.h"

#include <cassert>


// These classes are implemented afterwards for readability
class Fitness_Fluct;
class Fitness_Phenotype;
class Fitness_Stability;

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
        static double compute(const Phenotype& phenotype, const Population &);
        static Phenovec current_optimum();
        		
		static Phenovec expand_vec(const Phenovec&, unsigned int);

    protected:
        // constructors/destructor
        Fitness(const ParameterSet&);    
    
        static Fitness * instance;
        
        Fitness_Fluct * fitfluct;
        Fitness_Phenotype * fitphen;
        Fitness_Stability * fitstab;
};



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
		
		virtual Phenovec get_new_strength(const Phenovec & old_strength, unsigned int generation) = 0;
		virtual Phenovec get_new_optimum(const Phenovec & old_optimum, unsigned int generation) = 0;
};

class Fitness_Fluct_Nofluct: public Fitness_Fluct
{
	public:
		Fitness_Fluct_Nofluct(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);		
};

class Fitness_Fluct_States: public Fitness_Fluct
{ // Abstract class
	public:
		Fitness_Fluct_States(const ParameterSet &);
		virtual ~Fitness_Fluct_States() = 0;  
		
	protected:
		Phenovec strength_state1;
		Phenovec strength_state2;
		Phenovec optimum_state1;
		Phenovec optimum_state2;
		unsigned int period;		
};

class Fitness_Fluct_Flips: public Fitness_Fluct_States
{
	public:
		Fitness_Fluct_Flips(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);
};

class Fitness_Fluct_Sflips: public Fitness_Fluct_States
{
	public: 
		Fitness_Fluct_Sflips(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);		
};

class Fitness_Fluct_Smooth: public Fitness_Fluct_States
{
	public:
		Fitness_Fluct_Smooth(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);	
};

class Fitness_Fluct_Noise: public Fitness_Fluct
{ // Abstract class
	public:
		Fitness_Fluct_Noise(const ParameterSet &);
		virtual ~Fitness_Fluct_Noise() = 0;	
			
	protected:
		Phenovec strength_sd;
		Phenovec optimum_sd;
		unsigned int period;
};

class Fitness_Fluct_Whitenoise: public Fitness_Fluct_Noise
{
	public:
		Fitness_Fluct_Whitenoise(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);		
	protected:
		Phenovec strength_ref;
		Phenovec optimum_ref;
};

class Fitness_Fluct_Brownian: public Fitness_Fluct_Noise
{
	public:
		Fitness_Fluct_Brownian(const ParameterSet &);
		Phenovec get_new_strength(const Phenovec &, unsigned int);
		Phenovec get_new_optimum(const Phenovec &, unsigned int);
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
	
		double get_fitness(const Phenotype & phenotype, const Population &);
		virtual double get_fitness_trait(unsigned int trait, const Phenotype & phenotype, const Population&) = 0;
		virtual void update(const Population &) { }
		virtual void fluctuate(Fitness_Fluct *, unsigned int) { }
		virtual Phenovec get_optimum() const;
};

class Fitness_Phenotype_Noselection: public Fitness_Phenotype
{
	public:
		Fitness_Phenotype_Noselection(const ParameterSet &);
		double get_fitness_trait(unsigned int trait, const Phenotype&, const Population&) { return(1.0); }
};

class Fitness_Phenotype_Directional: public Fitness_Phenotype
{ // Abstract class for directional selection
	public: 
		Fitness_Phenotype_Directional(const ParameterSet &);
		virtual ~Fitness_Phenotype_Directional() = 0;
		virtual double get_fitness_trait(unsigned int trait, const Phenotype &, const Population &) = 0;
		void update(const Population &);
		void fluctuate(Fitness_Fluct *, unsigned int);
	protected:
		Phenovec strength;
		Phenovec popmean;
		const Population * popmem;
};

class Fitness_Phenotype_Linear: public Fitness_Phenotype_Directional
{
	public: 
		Fitness_Phenotype_Linear(const ParameterSet &);
		double get_fitness_trait(unsigned int trait, const Phenotype &, const Population&);
};

class Fitness_Phenotype_Expo: public Fitness_Phenotype_Directional
{
	public: 
		Fitness_Phenotype_Expo(const ParameterSet &);
		double get_fitness_trait(unsigned int trait, const Phenotype &, const Population&);
};

class Fitness_Phenotype_Concave: public Fitness_Phenotype_Directional
{
	public:
		Fitness_Phenotype_Concave(const ParameterSet &);
		double get_fitness_trait(unsigned int trait, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Stabilizing: public Fitness_Phenotype
{ // Abstract class for stabilizing selection
	public:
		Fitness_Phenotype_Stabilizing(const ParameterSet &);
		virtual ~Fitness_Phenotype_Stabilizing() { }
		virtual double get_fitness_trait(unsigned int, const Phenotype&, const Population&) = 0;
		void fluctuate(Fitness_Fluct *, unsigned int);
		Phenovec get_optimum() const { return(optimum);}
				
	protected:
		Phenovec strength;
		Phenovec optimum;
};

class Fitness_Phenotype_Gaussian: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Gaussian(const ParameterSet &);
		double get_fitness_trait(unsigned int, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Quadratic: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Quadratic(const ParameterSet &);
		double get_fitness_trait(unsigned int, const Phenotype&, const Population&);
};

class Fitness_Phenotype_Biconvex: public Fitness_Phenotype_Stabilizing
{
	public:
		Fitness_Phenotype_Biconvex(const ParameterSet &);
		double get_fitness_trait(unsigned int, const Phenotype&, const Population&);
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
		
		double get_fitness(const Phenotype &);
		virtual double get_fitness_trait(unsigned int, const Phenotype &) = 0;
};

class Fitness_Stability_Noselection: public Fitness_Stability
{
	public:
		Fitness_Stability_Noselection(const ParameterSet &);
		double get_fitness_trait(unsigned int, const Phenotype &);	
};

class Fitness_Stability_Exponential: public Fitness_Stability
{
	public:
		Fitness_Stability_Exponential(const ParameterSet &);
		double get_fitness_trait(unsigned int, const Phenotype &);
	protected:
		Phenovec strength;
};

#endif // FITNESS_H_INCLUDED

