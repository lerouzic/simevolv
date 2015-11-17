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
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir) const;
	    
	protected:
	    ArchiAdditive() { }	
	
	private:
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive &, const unsigned int);
};

template<class Archive>
void ArchiAdditive::serialize(Archive & ar, const unsigned int version)
{
	ar & boost::serialization::base_object<Architecture>(*this);
	// nothing to do here
}

#endif // ARCHIADDITIVE_H_INCLUDED
