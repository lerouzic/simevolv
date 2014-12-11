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


/* constructor called by Haplotype */
Allele::Allele(const vector<double> content)
	: allele(content)
{
}

Allele::Allele(const Allele & copy)
	: allele(copy.allele)
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


/* adds of the value of the two alleles of a genotype
 * (static function) */
vector<double> Allele::combine_add(const Allele & a1, const Allele & a2) 
{
	vector<double> ans;
	unsigned int all_size =  a1.allele.size();
	for (unsigned int sa = 0; sa < all_size; sa++) 
	{
		ans.push_back(a1.allele[sa] + a2.allele[sa]);
	}
	return(ans);
}

vector<double> Allele::combine_mean(const Allele & a1, const Allele & a2) 
{
	vector<double> ans;
	unsigned int all_size =  a1.allele.size();
	for (unsigned int sa = 0; sa < all_size; sa++) 
	{
		ans.push_back(0.5*(a1.allele[sa] + a2.allele[sa]));
	}
	return(ans);
}

shared_ptr<Allele> Allele::make_mutant(double mutsd) const
{
    int mutated_site = floor(allele.size()*Random::randnum());
    double modifier = mutsd * Random::randgauss();
    shared_ptr<Allele> a(new Allele(*this));
    a->allele[mutated_site] += modifier;
    return(a);
}



/******************** DERIVED CLASSES *********************************/
// Allele_zero

Allele_zero::Allele_zero(const vector<double> content)
	: Allele(content)
{
}

Allele_zero::Allele_zero(const Allele_zero & copy)
	: Allele(copy.allele)
{
}

shared_ptr<Allele> Allele_zero::make_mutant(double mutsd) const
{
	vector<unsigned int> non_zero;
	for (unsigned int i = 0; i < allele.size(); i++) {
		if (allele[i] != 0.0) 
			non_zero.push_back(i);
	}
	// This can't be a shared_ptr<Allele> since we need to access the protected member of the base
	// class. C++ is sometimes a bit too subtle... 
    shared_ptr<Allele_zero> a(new Allele_zero(*this));	 
	if (non_zero.size() > 0) {
		int mutated_site = non_zero[floor(non_zero.size()*Random::randnum())];
		double modifier = mutsd * Random::randgauss();
		a->allele[mutated_site] += modifier;		
	}
	// Note: if all sites are 0, no mutation, but
	// we create a new Allele instance.
	
    return(dynamic_pointer_cast<Allele>(a));
}
