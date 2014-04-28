// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Parconst.h"
#include "Allele.h"
#include "Architecture.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;



// constructors and destructor

/* default constructor */
Allele::Allele()
{
    int sall = Allele::all_size();
    
    cerr << "Calling Allele default constructor: should probably not happen." << endl;

    for(int i = 0; i < sall; i++)
    {
        allele.push_back(0.0);
    }
}


/* constructor called by Haplotype */
Allele::Allele(const vector<double> content)
	: allele(content)
{
}


// operator overload

int Allele::operator==(const Allele& other) const
{
    return((*this).allele == other.allele);
}


int Allele::operator!=(const Allele& other) const
{
    return(!(*this == other));
}


// functions

/* collect the size of the allele from the Architecture and return it */
unsigned int Allele::all_size() const
{
    Architecture * archi = Architecture::Get();
    int sall = archi -> all_size();
    return sall;
}


/* mean of the value of the two alleles of a genotype
 * (static function) */
vector<double> Allele::combine_add(const Allele & a1, const Allele & a2) 
{
	vector<double> ans;
	unsigned int all_size =  a1.allele.size();
	for (unsigned int sa = 0; sa < all_size; sa++) 
	{
		ans.push_back(0.5*(a1.allele[sa] + a2.allele[sa]));
	}
	return(ans);
}


