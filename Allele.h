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



#ifndef ALLELE_H_INCLUDED
#define ALLELE_H_INCLUDED

#include <vector>

#include "Parameters.h"


class Allele
{
    friend class Haplotype;
    friend class Architecture; 
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryWagner;
    friend class ArchiWagner;
    friend class ArchiMasel;
	
	public :
	    //constructors/destructor
	    Allele();
	    Allele(const std::vector<double>);
	
	    //operator overload
	    int operator== (const Allele&) const;
	    int operator!= (const Allele&) const;
	
	    //functions
	    unsigned int all_size() const;
		static std::vector<double> combine_add(const Allele &, const Allele &);
	
	protected :
	    std::vector<double> allele;
};


#endif // ALLELE_H_INCLUDED
