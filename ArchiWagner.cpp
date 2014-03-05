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



#include "ArchiWagner.h"
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <vector>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;


// constructors and destuctor

ArchiWagner::ArchiWagner()
{
    assert(false); // The default constructor should never be used.
}


ArchiWagner::ArchiWagner(const ParameterSet& param) 
    : ArchiRegulatoryWagner(param)
{
}


// functions

Phenotype ArchiWagner::phenotypic_value (const Genotype& genotype) const
{
	// creation of the W matrix;
	vector<vector<double> > matrix;
	for (unsigned int loc=0; loc < nb_loc(); loc++) 
	{
		matrix.push_back(Allele::combine_add(*genotype.gam_father.haplotype[loc], *genotype.gam_mother.haplotype[loc]));
	}
	
	// creation of the w_matrix and st_vector for using ublas 
	//using namespace boost::numeric::ublas;
	
	unsigned int nloc = nb_loc();
	boost::numeric::ublas::vector<double> St(nloc);
	boost::numeric::ublas::matrix<double> W(nloc, nloc); 
		
	for (unsigned int i=0 ; i<nloc ; i++)
	{
		St(i) = So[i];
	}
	
	for (unsigned int i=0; i<nloc; i++) 
	{
        for (unsigned int j=0; j<nloc; j++) 
        {
            W(i,j) = matrix[i][j];
        }
	}
		
	// simulation
	boost::numeric::ublas::vector<double> h(nloc);
	
	for (unsigned int t=0 ; t<timesteps ;t++)
	{
		h = prod(St,W);
		for (unsigned int i=0 ; i<h.size() ; i++)
		{
			if (h(i)<0) 
			{
				St(i) = -1.;
			}
			else if (h(i)>0) 
			{
				St(i) = 1.;
			}
			else 
			{
				St(i) = 0.;
			}
		}
	}
	
	// output
	std::vector<double> Sf;
	for (unsigned int i=0 ; i<nloc ; i++)
	{
		Sf.push_back(St(i));
	}
	
	return Phenotype(Sf);
}
