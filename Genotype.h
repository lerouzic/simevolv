// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

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
#include "Phenotype.h"

#include <iostream>


class Genotype
{
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiMasel;
    friend class ArchiSiegal; 
    friend class ArchiM2;
    friend class ArchiBoolean;

	public:
	    //constructors / destructors
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
	    Phenotype phenotypic_value(const Genotype&) const;
	    void draw_mutation();
	    void make_mutation(bool test = false);
	
	protected:
	    Haplotype gam_father;
	    Haplotype gam_mother;
};

#endif // GENOTYPE_H_INCLUDED
