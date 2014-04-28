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



#ifndef ARCHIREGULATORYWAGNER_H_INCLUDED
#define ARCHIREGULATORYWAGNER_H_INCLUDED

#include "Architecture.h"


class ArchiRegulatoryWagner : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiRegulatoryWagner();
	    ArchiRegulatoryWagner(const ParameterSet&);
	    ~ArchiRegulatoryWagner() {}
		
		// Inherited functions 
		virtual double sigma(double h) const {return(h)};
		virtual Phenotype phenotypic_value(const Genotype&) const;
		std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int) const;
	
	protected :
		unsigned int sall;
		std::vector<double> So;
		std::vector<std::vector<double> > connectivity_matrix; // this contains initial allelic values (for clonal pops), not only 0 or 1
		unsigned int timesteps;
		double basal;
		
	    //functions
		void init_connectivity_matrix(const ParameterSet &);		

};

#endif // ARCHIREGULATORYWAGNER_H_INCLUDED
