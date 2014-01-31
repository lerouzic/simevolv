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



#ifndef GENOTYPE_H_INCLUDED
#define GENOTYPE_H_INCLUDED

#include "Haplotype.h"
#include "Parameters.h"

#include <iostream>



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
	    Genotype(const ParameterSet &);
	
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


