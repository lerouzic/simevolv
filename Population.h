#ifndef POPULATION_H_INCLUDED
#define POPULATION_H_INCLUDED


#include "Parameters.h"
#include "Architecture.h"
#include "Individual.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>


class Individual;


class Population
{
public :
    //constructors/destructors
    Population();
    Population(long int);
    Population(const Population&);
    Population(const std::vector<Individual>&);
    Population(const ParameterSet&);

    //operator overload
    Population& operator = (const Population&);
    int operator == (const Population&) const;

    //instance/initialization
    void initialize(const ParameterSet &);

    //functions
    Population reproduce(long int offspr_number = 0) const;
    void update(void);
    //~ std::vector<double> phenotypes() const;
    double mean_phenotype() const;
    long int size() const;
    std::vector<double> cumul_fitness() const;
    const Individual & pick_parent(const std::vector<double>&) const;
    long int search_fit_table(double, const std::vector<double>&) const;
    long int sequential_search_fit_table(double, const std::vector<double>&) const;
    Individual iterator_search_fit_table(double, const std::vector<double>&) const;
    void draw_mutation();
    void make_mutation();

    //output
    void write() const;
    void write_debug(std::ostream&) const;
    void write_xml(std::ostream&) const;
    void write_simple(std::ostream&) const;
    void write_summary(std::ostream&) const;


protected :
    std::vector<Individual> pop;

};

#endif // POPULATION_H_INCLUDED
