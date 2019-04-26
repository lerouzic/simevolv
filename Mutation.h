// Copyright 2019 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef MUTATION_H_INCLUDED
#define MUTATION_H_INCLUDED

#include "types.h"
#include "Parameters.h"

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#endif


//////////////////////////////////////////////////////////////
// Mutation types

class MutType
{ // This is an interface
	public:
		MutType() { }
		// MutType(const MutType &) = delete;
		virtual ~MutType() { }
		
		virtual allele_type mutate(allele_type oldv) const = 0;
		virtual MutType * clone() const = 0;
		
	private:
		#ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) { }
		#endif
};

class MutBoolean: public MutType
{
	public:
		MutBoolean ();
		MutBoolean(const MutBoolean &);
		~MutBoolean() { }
	
		allele_type mutate(allele_type oldv) const;
		virtual MutType * clone() const;
		
	private:
		#ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<MutType>(*this);
		}
		#endif
};

class MutGaussianCumul: public MutType
{
	public:
		MutGaussianCumul(allele_type sd);
		MutGaussianCumul(const MutGaussianCumul &);
		~MutGaussianCumul() { }
		
		allele_type mutate(allele_type oldv) const;
		virtual MutType * clone() const;
		
	protected:
		allele_type mutsd;
		
	private:
		MutGaussianCumul() : mutsd(0.0) { }
		
		#ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<MutType>(*this);
			ar & mutsd;
		}
		#endif
};

class MutGaussianStationary: public MutType
{
	public:
		MutGaussianStationary(allele_type sd);
		MutGaussianStationary(const MutGaussianStationary &);
		~MutGaussianStationary() { }
		
		allele_type mutate(allele_type oldv) const;
		virtual MutType * clone() const;
	
	protected:
		allele_type mutsd;
		
	private:
		MutGaussianStationary() : mutsd(0.0) { }
	
		#ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{
			ar & boost::serialization::base_object<MutType>(*this);
			ar & mutsd;
		}
		#endif
};

///////////////////////////////////////////////
// Mutation model (the interface that should be used)


class MutationModel
{
	public:
		MutationModel(); // not really useful, but necessary for serialization. Use with caution. 
		MutationModel(const ParameterSet& param, unsigned int locus=0, bool test=false);
		MutationModel(const MutationModel&);
		~ MutationModel();
		MutationModel & operator= (const MutationModel &);
		
		
		allele_type mutate(allele_type oldv) const; 
		
	protected:
		MutType * mut;
		
	private:
		#ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version)
		{ // Why doesn't BOOST_CLASS_EXPORT do the job? 
			ar.template register_type<MutBoolean>();
			ar.template register_type<MutGaussianCumul>();
			ar.template register_type<MutGaussianStationary>();
			ar & mut;
		}
		#endif
};




#endif // MUTATION_H_INCLUDED
