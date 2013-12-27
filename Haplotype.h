#ifndef HAPLOTYPE_H_INCLUDED
#define HAPLOTYPE_H_INCLUDED


#include "Parameters.h"
#include "Allele.h"

#include <iostream>
#include <string>
#include <vector>
#include <numeric>



class Haplotype
{
    friend class Genotype;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;

public :
    //constructors/destructor
    Haplotype();

    //operator overload
    int operator== (const Haplotype&) const;
    int operator!= (const Haplotype&) const;

    //functions
    int nb_loc() const;
    int all_size() const;
    void print() const;
    void draw_mutation();
    void make_mutation();

    //output/debug
    void write_debug (std::ostream&) const;
    void write_xml   (std::ostream&) const;
    void write_simple(std::ostream&) const;

protected :
    std::vector<Allele> haplotype;
};


#endif // HAPLOTYPE_H_INCLUDED

