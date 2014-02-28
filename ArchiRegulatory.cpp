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



#include "ArchiRegulatory.h"
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

using namespace std;


// constructors and destuctor

ArchiRegulatory::ArchiRegulatory()
{
    assert(false); // The default constructor should never be used.
}

ArchiRegulatory::ArchiRegulatory(const ParameterSet& param) 
    : Architecture(param)
    , sall(nb_loc())
    , timesteps(4)
    , basal(0.1)
{
	/* Update when the parameter file will be OK */
	
	init_connectivity_matrix(param); // creates connectivity_matrix

	for (unsigned int n=0; n<nb_loc(); n++)
	{
		so.push_back(basal);
	}	
}

shared_ptr<Allele> ArchiRegulatory::allele_init(const ParameterSet & param, unsigned int loc) const
{
	bool clonal = (param.getpar(INIT_CLONAL)->GetString() == CL_clonal);
	vector<double> temp_allele;
	for (unsigned int i = 0; i < sall; i++) {
		assert(connectivity_matrix.size() > loc);
		assert(connectivity_matrix[loc].size() > i);
		if (connectivity_matrix[loc][i] == 0.0) {
			temp_allele.push_back(0.0);
		} else {
			double value;
			if (clonal) {
				value = connectivity_matrix[loc][i];
			} else {
				value = param.getpar(INIT_ALLELES) -> GetDouble();
			}
			temp_allele.push_back(value);
		}
	}
	shared_ptr<Allele> a (new Allele(temp_allele));
	return(a);
}


// (protected) functions

static unsigned int count_init = 0; // debug instruction

void ArchiRegulatory::init_connectivity_matrix(const ParameterSet & param)
{
	count_init++; // debug instruction;
	
	double connectivity = 0.5;
	//  = param.getpar(INIT_CONNECT)->GetDouble();
	// bool clonal = (param.getpar(INIT_CLONAL)->GetString() == CL_clonal);
	
	for (unsigned int loc = 0; loc < nb_loc(); loc++) {
		vector<double> allele_pattern;
		for (unsigned int n=0 ; n<sall ; n++)
		{
			if (Random::randnum() < connectivity)
			{
				// This has some interest only in clonal populations.
				// In non-clonal pops, this value will be overwritten anyway.
				// double value = 1.0 // in non-clonal, this would be enough.
				double value = param.getpar(INIT_ALLELES) -> GetDouble();
				allele_pattern.push_back(value);

			}
			else
			{
				allele_pattern.push_back(0);
			}
		}
		connectivity_matrix.push_back(allele_pattern);		
	}
	std::cout << count_init << "\n";
}



Phenotype ArchiRegulatory::phenotypic_value (const Genotype& genotype) const
{
	//~ ParameterSet param;
	//~ 
	//~ vector<vector<double> > matrix;
	//~ for (int i=0 ; i<nb_loc() ; i++)
	//~ {
		//~ matrix.push_back(init_pattern());
	//~ }
	//~ 
	//~ vector<shared_ptr<Allele> > pmatrix;
	//~ for (int i=0 ; i<nb_loc() ; i++)
	//~ {
		//~ pmatrix.push_back(Allele(init_pattern()));
	//~ }
//~ 
//~ 
	//~ Haplotype haplotype;
	//~ for (int i=0 ; i<(param.getpar(INIT_PSIZE)->GetInt()) ; i++)
	//~ {
		//~ for (int g=0 ; g<2 ; g++)
		//~ {
			//~ for (int n=0 ; n<nb_loc() ; n++)
			//~ {
				//~ haplotype.push_back(pmatrix[n]);
			//~ } 
		//~ }
	//~ }
	//~ 
	//~ return (haplotype[1]);
	//~ 
	//~ 
	//~ 
	//~ 
	//~ 
	//~ 
	//~ vector<int> st;
	//~ vector<int> h;
	//~ for (int t=0 ; t<timesteps ;t++)
	//~ {
		//~ if (t==0)
		//~ {
			//~ h = so * w;
			//~ if (h<0) {st = -1 * h;}
			//~ else if (h>0) {st = 1 * h;}
			//~ else {st = O * h;}
		//~ }
		//~ else
		//~ {
			//~ h = st * w;
			//~ if (h<0) {st = -1 * h;}
			//~ else if (h>0) {st = 1 * h;}
			//~ else {st = O * h;}
		//~ }
	//~ }
	//~ return st;

}
