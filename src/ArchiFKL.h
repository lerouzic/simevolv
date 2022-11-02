// Copyright 2022       Arnaud Le Rouzic    <arnaud.le-rouzic@universite-paris-saclay.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


#ifndef ARCHIFKL_H_INCLUDED
#define ARCHIFKL_H_INCLUDED

#include "Architecture.h"
#include "ArchiAdditive.h"
#include "Allele.h"

#include <iostream>


class ArchiFKL : public ArchiAdditive
{
	public :
		//constructors/destructor
		ArchiFKL(const Architecture&) = delete;
		ArchiFKL(const ParameterSet&);
		virtual ~ArchiFKL() { }
		
		virtual unsigned int nb_phen() const; 
		
		virtual std::vector<rate_type> mutation_rates(const Haplotype &) const;
		virtual std::vector<rate_type> mutmutation_rates() const;
	
		virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo &, bool sdinittest = false, bool sddynamtest = false) const;
		virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int loc = 0) const;
		virtual std::shared_ptr<Allele> allele_mutation(std::shared_ptr<Allele>, unsigned int loc = 0, bool test = false) const;
		virtual std::shared_ptr<Allele> allele_mut_mutation(std::shared_ptr<Allele>, unsigned int loc = 0) const;

	protected :
		unsigned int nphen;
		std::vector<rate_type> mutmutrate;

		ArchiFKL() { }; // Necessary for serialization
		virtual void update_param_internal(const ParameterSet&);
		
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
			ar & boost::serialization::base_object<ArchiAdditive>(*this);
			ar & mutmutrate;
		}
        #endif
};


#endif // ARCHIFKL_H_INCLUDED
