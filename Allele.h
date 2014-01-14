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



class Allele
{
    friend class Haplotype;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
	
	public :
	    //constructors/destructor
	    Allele(int nall);
	
	    //operator overload
	    int operator== (const Allele&) const;
	    int operator!= (const Allele&) const;
	
	    //functions
	    int all_size() const;
	    void print() const;
	    void make_mutation(int);
	
	protected :
	    std::vector<double> allele;
};


#endif // ALLELE_H_INCLUDED
