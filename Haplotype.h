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

#include <iostream>
#include <vector>



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

