// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef HAPLOTYPE_H_INCLUDED
#define HAPLOTYPE_H_INCLUDED

#include "Allele.h"
#include "Parameters.h"

#include <iostream>
#include <vector>
#include <memory> // shared_ptr



class Haplotype
{
    friend class Genotype;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryWagner;
    friend class ArchiWagner;
    friend class ArchiMasel;  

public :
    //constructors/destructor
    Haplotype();
    Haplotype(const ParameterSet &);
    Haplotype(const Haplotype &);
    Haplotype(const std::vector<std::shared_ptr<Allele> > &);

    //operator overload
    int operator== (const Haplotype&) const;
    int operator!= (const Haplotype&) const;

    //functions
    unsigned int nb_loc() const;
    void draw_mutation();				// draws a mutation according to the mutation rate
    void make_mutation(); 				// makes a mutation at a random locus
    void make_mutation(unsigned int); 	// makes a mutation at a specific locus

    //output/debug
    void write_debug (std::ostream&) const;
    void write_xml   (std::ostream&) const;
    void write_simple(std::ostream&) const;

protected :
    std::vector<std::shared_ptr<Allele> > haplotype;
};


#endif // HAPLOTYPE_H_INCLUDED

