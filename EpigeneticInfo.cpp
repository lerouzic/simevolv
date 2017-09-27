// Copyright 2007-2016 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>
// Copyright 2016      Andreas Odorico     

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "EpigeneticInfo.h"
#include "Individual.h"

#include <cassert>

using namespace std;

EpigeneticInfo::EpigeneticInfo() 
	: epigenet(0.0)
	// This constructor exists for convenience, but the object
	// does not contain any information. 
{ }

EpigeneticInfo::EpigeneticInfo(const EpigeneticInfo & copy)
	: epigenet(copy.epigenet)
	, mother_phen(copy.mother_phen)
{
	assert(epigenet >= 0.0);
	assert(epigenet <= 1.0);	
}

EpigeneticInfo::EpigeneticInfo(const double epi, const Phenotype & mothphen)
	: epigenet(epi)
	, mother_phen(mothphen)
{ 
	assert(epigenet >= 0.0);
	assert(epigenet <= 1.0);
}

EpigeneticInfo::EpigeneticInfo(const Individual & mother)
	: epigenet(mother.get_epigenet())
	, mother_phen(mother.get_phenotype())
{ 
	assert(epigenet >= 0.0);
	assert(epigenet <= 1.0);	
}

EpigeneticInfo::~EpigeneticInfo() 
{ }

double EpigeneticInfo::get_epigenet() const
{
	return(epigenet);
}

EpigeneticInfo & EpigeneticInfo::operator = (const EpigeneticInfo& copy) 
{
	if (this == &copy)
		return (*this);
		
	epigenet = copy.epigenet;
	mother_phen = copy.mother_phen;
	
	return(*this);
}

Phenotype EpigeneticInfo::get_phenotype() const
{
	return(mother_phen);
}

bool EpigeneticInfo::is_defined() const
{
	return(mother_phen.is_defined());
}

