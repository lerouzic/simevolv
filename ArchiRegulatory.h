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
	    ArchiRegulatory(const ParameterSet&);
	    ~ArchiRegulatory() {}
		
		// Inherited functions
	    Phenotype phenotypic_value(const Genotype&) const;
		std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
	
	protected :
		unsigned int sall;
		std::vector<double> so;
		std::vector<std::vector<double> > connectivity_matrix; // this contains initial allelic values (for clonal pops), not only 0 or 1
		unsigned int timesteps;
		double basal;
		
	    //functions
		void init_connectivity_matrix(const ParameterSet &);		

};

#endif // ARCHIREGULATORY_H_INCLUDED
