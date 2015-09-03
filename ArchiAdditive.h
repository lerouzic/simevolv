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
	    ArchiAdditive() = delete;
	    ArchiAdditive(const Architecture&) = delete;
	    ArchiAdditive(const ParameterSet&);
	    virtual ~ArchiAdditive() {}
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    //functions
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir) const;
};

#endif // ARCHIADDITIVE_H_INCLUDED
