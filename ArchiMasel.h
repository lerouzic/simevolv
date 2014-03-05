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



#ifndef ARCHIMASEL_H_INCLUDED
#define ARCHIMASEL_H_INCLUDED

#include "ArchiRegulatoryWagner.h"



class ArchiMasel : public ArchiRegulatoryWagner
{
	public :
	    //constructors/destructor
	    ArchiMasel();
	    ArchiMasel(const ParameterSet&);
	    ~ArchiMasel() {}
		
		// Inherited functions
	    Phenotype phenotypic_value(const Genotype&) const;
			
	protected :
		
};

#endif // ARCHIMASEL_H_INCLUDED
