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



#ifndef ARCHIADDITIVE_H_INCLUDED
#define ARCHIADDITIVE_H_INCLUDED

#include "Architecture.h"

#include <iostream>

class ArchiAdditive : public Architecture
{	
	public :
	    //constructors/destructor
	    ArchiAdditive(const Architecture&) = delete;
	    ArchiAdditive(const ParameterSet&);
	    virtual ~ArchiAdditive();
	
	    //functions
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir, const EpigeneticInfo &, bool sdinittest = false, bool sddynamtest = false) const;
		// virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0, bool test = false) const; 
		// the default is acceptable

	protected:
	    ArchiAdditive() { }	
	
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & boost::serialization::base_object<Architecture>(*this);
        }
        #endif
};


#endif // ARCHIADDITIVE_H_INCLUDED
