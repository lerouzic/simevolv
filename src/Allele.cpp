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



#include "Allele.h"

#include "Parconst.h"
#include "Architecture.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <cmath>

 #ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT(Allele)
#endif

using namespace std;

// constructors and destructor

Allele::Allele(const vector<allele_type> content)
	: allele(content)
	, type_allele(vector<string>(content.size(), TA_norm))
{
}

Allele::Allele(const vector<allele_type> content, const vector<string> all_type)
	: allele(content)
	, type_allele(all_type)
{
}

Allele::Allele(const Allele & copy)
	: allele(copy.allele)
	, type_allele(copy.type_allele)
{
}

// operator overload
bool Allele::operator==(const Allele& other) const
{
    return((*this).allele == other.allele);
}

bool Allele::operator!=(const Allele& other) const
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
vector<allele_type> Allele::combine_add(const Allele & a1, const Allele & a2) 
{
	vector<allele_type> ans;
	unsigned int all_size =  a1.allele.size();
	for (unsigned int sa = 0; sa < all_size; sa++) 
	{
		ans.push_back(a1.allele[sa] + a2.allele[sa]);
	}
	return(ans);
}

vector<allele_type> Allele::combine_mean(const Allele & a1, const Allele & a2) 
{
	vector<allele_type> ans;
	unsigned int all_size =  a1.allele.size();
	for (unsigned int sa = 0; sa < all_size; sa++) 
	{
		ans.push_back(0.5*(a1.allele[sa] + a2.allele[sa]));
	}
	return(ans);
}

//Allele combination used in the Boolean Architecture. Combines the two alleles by using an OR function
vector<allele_type> Allele::combine_OR(const Allele & a1, const Allele & a2)
{
    vector<allele_type> ans;
    
    unsigned int all_size =  a1.allele.size();
    for (unsigned int sa = 0; sa < all_size; sa++)
    {
        if(a1.allele[sa]||a2.allele[sa])
            ans.push_back(1);
        else
            ans.push_back(0);
    }
    return(ans);
}

shared_ptr<Allele> Allele::make_mutant_at_site(size_t site, const MutationModel& mutmodel) const
{
	// Creates a new allele even if the site is not mutable. To be avoided. 
	shared_ptr<Allele> a(new Allele(*this));
	a->allele[site] = mutmodel.mutate(a->allele[site], a->type_allele[site]);
	return(a);
}

shared_ptr<Allele> Allele::make_mutant_random_site(const MutationModel& mutmodel) const
{
	// Here the algorithm is more complex, as we don't want to risk a mutation in a non-mutable site
	vector<size_t> non_zero;
	for (size_t i = 0; i < allele.size(); i++) {
		if ( !(( type_allele[i] == TA_immut) || ((type_allele[i] == TA_zero && allele[i] == 0.0))) ) 
			non_zero.push_back(i);
	}
	if (non_zero.size() == 0) {
		// Not ideal, a new shared_ptr is created while its content is identical to *this. 
		return(make_shared<Allele>(*this)); 
	}
	size_t mutated_site = non_zero.at(floor(non_zero.size()*Random::randnum()));
	return(make_mutant_at_site(mutated_site, mutmodel));
}

shared_ptr<Allele> Allele::make_mutant_all_sites(const MutationModel& mutmodel) const
{
	shared_ptr<Allele> a(new Allele(*this));
	for(size_t i = 0; i < allele.size(); i++) {
		a->allele[i] = mutmodel.mutate(a->allele[i], a->type_allele[i]);
	}
	return(a);
}


std::vector<allele_type> Allele::get_raw() const
{
	return(allele);
}
