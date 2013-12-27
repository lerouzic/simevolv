#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include "Parameters.h"
#include "Architecture.h"
#include "Genotype.h"
#include "Population.h"

#include <iostream>
#include <vector>



class Population;


class Individual
{
public :
    //constructors/destructor
    Individual();
    Individual(const Haplotype&, const Haplotype&);
    Individual(const Individual&);
    Individual(const ParameterSet&);
    virtual ~Individual();

    // operator overload
    Individual & operator = (const Individual&);
    int operator == (const Individual&) const;

    // instance/initialization
    void initialize();

    //functions
    void update_fitness(const Population &);
    void update_fitness(const double);
    double get_fitness() const;
    double get_phenotype() const;
    Haplotype produce_gamete() const;
    static Individual mate(const Individual&, const Individual&);
    void draw_mutation();
    void make_mutation();

    //output
    void write_debug (std::ostream&) const;
    void write_xml (std::ostream&) const;
    void write_simple(std::ostream&) const;


public :
    Genotype genotype;
    double phenotype;
    double fitness;

};


#endif // INDIVIDUAL_H_INCLUDED
