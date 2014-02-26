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



#ifndef ARCHIREGULATORY_H_INCLUDED
#define ARCHIREGULATORY_H_INCLUDED

#include "Architecture.h"



class ArchiRegulatory : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatory();
	    ArchiRegulatory(const Architecture&);
	    ArchiRegulatory(const ParameterSet&);
	    ~ArchiRegulatory() {}
		
		//operator overload

	    //functions
	    double init_value() const;
	    std::vector<double> create_pattern() const;
	    Phenotype phenotypic_value(const Genotype&) const;
	
	protected :
		int sall;
		int init;
		std::vector<double> so;
		int timesteps;
		double basal;
		double connectivity;

};

#endif // ARCHIREGULATORY_H_INCLUDED
