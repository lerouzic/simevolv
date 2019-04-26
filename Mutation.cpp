// Copyright 2019 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>


/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Mutation.h"

#include "Parconst.h"
#include "Random.h"

#include <cassert>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>

BOOST_SERIALIZATION_ASSUME_ABSTRACT(MutType);
BOOST_CLASS_EXPORT(MutBoolean);
BOOST_CLASS_EXPORT(MutGaussianCumul);
BOOST_CLASS_EXPORT(MutGaussianStationary);
#endif


MutationModel::MutationModel()
	: mut(nullptr)
{ // Nothing can be done whith the default-constructed object, but this is necessary to store MutationModels in a vector
}

MutationModel::MutationModel(const MutationModel & mm)
	: mut (mm.mut->clone())
{
}

MutationModel::MutationModel(const ParameterSet & param, unsigned int locus /* = 0 */, bool test /* = false */)
	: mut(nullptr)
{
	// The task of the constructor is to initialize the mut pointer
	// to the right MutType class
	
	if (param.getpar(TYPE_ARCHI)->GetString() == AR_Boolean) {
		mut = new MutBoolean();
	} else {
		auto psd = test ? param.getpar(OUT_CANAL_MUTSD) : param.getpar(GENET_MUTSD);

		if (param.getpar(GENET_MUTMEM)->GetString() == MM_cumul) {
			mut = new MutGaussianCumul(psd->GetDouble(locus));
		} else if (param.getpar(GENET_MUTMEM)->GetString() == MM_stationary) {
			mut = new MutGaussianStationary(psd->GetDouble(locus));
		} else {
			std::cerr << "Impossible to initialize mutations" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	assert(mut);
}

MutationModel::~MutationModel()
{
	if (mut) {
		delete mut;
		mut = nullptr;
	}
}

MutationModel & MutationModel::operator= (const MutationModel & mm)
{
	if (mut) 
		delete mut;
	mut = nullptr;
	if (mm.mut)
		mut = mm.mut->clone();
	return *this;
}

allele_type MutationModel::mutate(allele_type oldv) const
{
	assert(mut);
	return mut->mutate(oldv); 
}



/////////////////// Mutation types /////////////////////

MutBoolean::MutBoolean()
{ // nothing to do here
}

MutBoolean::MutBoolean(const MutBoolean & mm)
{
}

allele_type MutBoolean::mutate(allele_type oldv) const
{
	assert(oldv >= 0. && oldv <= 1.);
	return 1. - oldv;
}

MutType * MutBoolean::clone() const
{
	return new MutBoolean(*this);
}



MutGaussianCumul::MutGaussianCumul(allele_type sd)
	: mutsd(sd)
{
	assert(mutsd >= 0.);
}

MutGaussianCumul::MutGaussianCumul(const MutGaussianCumul & mm)
	: mutsd(mm.mutsd)
{
}

allele_type MutGaussianCumul::mutate(allele_type oldv) const
{
	return oldv + mutsd * Random::randgauss();
}

MutType * MutGaussianCumul::clone() const
{
	return new MutGaussianCumul(*this);
}



MutGaussianStationary::MutGaussianStationary(allele_type sd)
	: mutsd(sd)
{
	assert(mutsd >= 0.);
}

MutGaussianStationary::MutGaussianStationary(const MutGaussianStationary & mm)
	: mutsd(mm.mutsd)
{
}

allele_type MutGaussianStationary::mutate(allele_type oldv) const
{
	return mutsd * Random::randgauss();
}

MutType * MutGaussianStationary::clone() const
{
	return new MutGaussianStationary(*this);
}

