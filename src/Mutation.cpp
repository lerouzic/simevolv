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
BOOST_CLASS_EXPORT(MutMultivGaussianCumul);
BOOST_CLASS_EXPORT(MutMultivGaussianStationary);

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
	
	std::string type_archi   = param.getpar(TYPE_ARCHI)  ->GetString();
	std::string genet_mutmem = param.getpar(GENET_MUTMEM)->GetString();
	
	if (type_archi == AR_Boolean) {
		mut = std::unique_ptr<MutType>(new MutBoolean());
	} else {
		auto psd = test ? param.getpar(OUT_CANAL_MUTSD) : param.getpar(GENET_MUTSD);
		
		if (type_archi == AR_FKL) {
			unsigned int numtrait = param.getpar(GENET_NBPHEN)->GetInt();
			
			std::vector<allele_type> mutsd;
			for (size_t i = 0; i < numtrait; i++)
				mutsd.push_back(psd->GetDouble(locus*numtrait + i));
				
			std::vector<allele_type> mutcor;
			for (size_t i = 0; i < numtrait - 1; i++)
				for (size_t j = i+1; j < numtrait; j++)
					mutcor.push_back(param.getpar(GENET_MUTCOR)->GetDouble(locus*numtrait*(numtrait-1)/2 + i*numtrait - i*(i+1)/2 + (j-i-1)));
			
			if (genet_mutmem == MM_cumul) {
				mut = std::unique_ptr<MutType> (new MutMultivGaussianCumul(mutsd, mutcor, param.getpar(GENET_MUTMUTSD)->GetDouble(locus)));
			} else if (genet_mutmem == MM_stationary) {
				mut = std::unique_ptr<MutType> (new MutMultivGaussianStationary(mutsd, mutcor, param.getpar(GENET_MUTMUTSD)->GetDouble(locus)));
			} else {
				std::cerr << "Impossible to initialize mutations" << std::endl;
				exit(EXIT_FAILURE);
			}
		} else {
			if (genet_mutmem == MM_cumul) {
				mut = std::unique_ptr<MutType> (new MutGaussianCumul(psd->GetDouble(locus)));
			} else if (genet_mutmem == MM_stationary) {
				mut = std::unique_ptr<MutType> (new MutGaussianStationary(psd->GetDouble(locus)));
			} else {
				std::cerr << "Impossible to initialize mutations" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	assert(mut);
}

MutationModel::~MutationModel()
{
}

MutationModel & MutationModel::operator= (const MutationModel & mm)
{
	mut = mm.mut->clone();
	return *this;
}

allele_type MutationModel::mutate(allele_type oldv, std::string type_all, size_t site /* = 0 */) const
{
	assert(mut);
	allele_type candidate = mut->mutate(oldv, site);
	if ((type_all == TA_immut) || ((type_all == TA_zero) && (oldv == 0.0)))
		return(oldv);
	if (type_all == TA_sign)
		return((((oldv > 0) && (candidate > 0)) || ((oldv < 0) && (candidate < 0))) ? candidate : -candidate);
	return(candidate);
}

std::vector<allele_type> MutationModel::mutate(std::vector<allele_type> oldv, std::vector<std::string> type_all) const
{
	assert(mut);
	
	std::vector<allele_type> candidate = mut->mutate(oldv);
	
	for (size_t i = 0; i < candidate.size(); i++) {
		if ((type_all[i] == TA_immut) || ((type_all[i] == TA_zero) && (oldv[i] == 0.0)))
			candidate[i] = oldv[i];
		if ((type_all[i] == TA_sign) && (((oldv[i] > 0) && candidate[i] < 0) || ((oldv[i] < 0) && (candidate[i] > 0)) ))
			candidate[i] = -candidate[i];
	}
	return(candidate);
}


rate_type MutationModel::mutate_mut(rate_type oldv) const
{
	assert(mut);
	auto candidate = mut->mutate_mut(oldv);
	return(candidate);
}

/////////////////// Mutation types /////////////////////

MutType::MutType()
{
}

MutType::~MutType()
{
}

std::vector<allele_type> MutType::mutate(std::vector<allele_type> oldv) const 
{ // If not specified in the derived function: just call the univariate version several times
	for (auto &oo : oldv)
		oo = this->mutate(oo);
	return oldv;
}






MutBoolean::MutBoolean()
{ // nothing to do here
}

MutBoolean::MutBoolean(const MutBoolean & mm)
{
}

allele_type MutBoolean::mutate(allele_type oldv, size_t site = 0) const
{
	assert(oldv >= 0. && oldv <= 1.);
	return 1. - oldv;
}

std::unique_ptr<MutType> MutBoolean::clone() const
{
	return std::unique_ptr<MutType> (new MutBoolean(*this));
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

allele_type MutGaussianCumul::mutate(allele_type oldv, size_t site = 0) const
{
	return oldv + mutsd * Random::randgauss();
}

std::unique_ptr<MutType> MutGaussianCumul::clone() const
{
	return std::unique_ptr<MutType> (new MutGaussianCumul(*this));
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

allele_type MutGaussianStationary::mutate(allele_type oldv, size_t site = 0) const
{
	return mutsd * Random::randgauss();
}

std::unique_ptr<MutType> MutGaussianStationary::clone() const
{
	return std::unique_ptr<MutType> (new MutGaussianStationary(*this));
}




MutMultivGaussianCumul::MutMultivGaussianCumul(std::vector<allele_type> sd, std::vector<allele_type> cor, allele_type mmutsd)
	: mutsd(sd)
	, mutcor(cor)
	, mutmutsd(mmutsd)
{
	// Temporary (?) limitation? So far, only bivariate distributions are handled
	assert(mutsd.size()  == 2);
	assert(mutcor.size() == 1);
	for (auto s : mutsd)
		assert(s >= 0.0);
	assert(mutmutsd > 0.0);
	for (auto r : mutcor)
		assert((r >= -1.0) && (r <= 1.0));
}

MutMultivGaussianCumul::MutMultivGaussianCumul(const MutMultivGaussianCumul & mm)
	: mutsd (mm.mutsd)
	, mutcor(mm.mutcor)
	, mutmutsd(mm.mutmutsd)
{
}

allele_type MutMultivGaussianCumul::mutate(allele_type oldv, size_t site = 0) const
{ // Does it really make sense? There will not be any change in the correlated sites
	return oldv + mutsd[site] * Random::randgauss();
}

std::vector<allele_type> MutMultivGaussianCumul::mutate(std::vector<allele_type> oldv) const
{
	std::vector<allele_type> multirand = Random::randbivgauss(mutcor[0]); // only works with bivariate distributions
	assert(oldv.size() == multirand.size());
	for (size_t i = 0; i < multirand.size(); i++)
		multirand[i] = oldv[i] + mutsd[i] * multirand[i];
	return multirand;
}

rate_type MutMultivGaussianCumul::mutate_mut(rate_type oldv) const
{
	return oldv + mutmutsd * Random::randgauss();
}

std::unique_ptr<MutType> MutMultivGaussianCumul::clone() const
{
	return std::unique_ptr<MutType> (new MutMultivGaussianCumul(*this));
}




MutMultivGaussianStationary::MutMultivGaussianStationary(std::vector<allele_type> sd, std::vector<allele_type> cor, allele_type mmutsd)
	: mutsd(sd)
	, mutcor(cor)
	, mutmutsd(mmutsd)
{
	// Temporary (?) limitation? So far, only bivariate distributions are handled
	assert(mutsd.size()  == 2);
	assert(mutcor.size() == 1);
	for (auto s : mutsd)
		assert(s >= 0.0);
	assert(mutmutsd >= 0.0);
	for (auto r : mutcor)
		assert((r >= -1.0) && (r <= 1.0));
}

MutMultivGaussianStationary::MutMultivGaussianStationary(const MutMultivGaussianStationary & mm)
	: mutsd (mm.mutsd)
	, mutcor(mm.mutcor)
	, mutmutsd(mm.mutmutsd)
{
}

allele_type MutMultivGaussianStationary::mutate(allele_type oldv, size_t site = 0) const
{
	return mutsd[site] * Random::randgauss();
}

std::vector<allele_type> MutMultivGaussianStationary::mutate(std::vector<allele_type> oldv) const
{ // no need for oldv, just here for the interface

	std::vector<allele_type> multirand = Random::randbivgauss(mutcor[0]); // only works with bivariate distributions

	for (size_t i = 0; i < multirand.size(); i++)
		multirand[i] = mutsd[i] * multirand[i];
	return multirand;
}

rate_type MutMultivGaussianStationary::mutate_mut(rate_type oldv) const
{
	return mutmutsd * Random::randgauss();
}


std::unique_ptr<MutType> MutMultivGaussianStationary::clone() const
{
	return std::unique_ptr<MutType> (new MutMultivGaussianStationary(*this));
}
