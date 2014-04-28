// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
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



class ArchiAdditive : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiAdditive();
	    ArchiAdditive(const Architecture&);
	    ArchiAdditive(const ParameterSet&);
	    ~ArchiAdditive() {}
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    //functions
	    Phenotype phenotypic_value(const Genotype&) const;
};

#endif // ARCHIADDITIVE_H_INCLUDED
