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

#ifndef EPIGENETICINFO_H_INCLUDED
#define EPIGENETICINFO_H_INCLUDED

#include "Phenotype.h"

#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#endif


class Individual;

class EpigeneticInfo {
	public:
		EpigeneticInfo();
		EpigeneticInfo(const EpigeneticInfo &);
		EpigeneticInfo(const double, const Phenotype &);
		EpigeneticInfo(const Individual &);
		~ EpigeneticInfo();
		EpigeneticInfo & operator = (const EpigeneticInfo &);
			
		double get_epigenet() const;
		Phenotype get_phenotype() const;
		bool is_defined() const;
	
	protected:
		double epigenet;
		Phenotype mother_phen;
        
	private:
    #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & epigenet;
            ar & mother_phen;
        }          
    #endif
};

#endif // EPIGENETICINFO_H_INCLUDED
