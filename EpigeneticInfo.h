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

class Individual;

class EpigeneticInfo {
	public:
		EpigeneticInfo();
		EpigeneticInfo(const EpigeneticInfo &);
		EpigeneticInfo(const double, const Phenotype &);
		EpigeneticInfo(const Individual &);
		virtual ~ EpigeneticInfo();
		EpigeneticInfo & operator = (const EpigeneticInfo) = delete;
			
		virtual double get_epigenet() const;
		virtual Phenotype get_phenotype() const;
		virtual bool is_defined() const;
	
	protected:
		const double epigenet;
		const Phenotype mother_phen;
		
};

#endif // EPIGENETICINFO_H_INCLUDED
