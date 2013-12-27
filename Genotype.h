#ifndef GENOTYPE_H_INCLUDED
#define GENOTYPE_H_INCLUDED


#include "Parameters.h"
#include "Haplotype.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>



class Genotype
{
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;

public:
    //constructors / destructors
    Genotype();
    Genotype(const Haplotype&, const Haplotype&);
    Genotype(const Genotype&);

    //operator overload
    int operator == (const Genotype&) const;
    int operator != (const Genotype&) const;

    //functions
    int nb_loc() const;
    int all_size() const;
    Haplotype recombine() const;
    double phenotypic_value(const Genotype&) const;
    void draw_mutation();
    void make_mutation();

    //output
    void write_debug (std::ostream&) const;
    void write_xml (std::ostream&) const;
    void write_simple(std::ostream&) const;

protected:
    Haplotype gam_father;
    Haplotype gam_mother;
};


#endif // GENOTYPE_H_INCLUDED


