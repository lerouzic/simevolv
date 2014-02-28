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


ArchiRegulatory::ArchiRegulatory(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


ArchiRegulatory::ArchiRegulatory(const ParameterSet& param) 
    : Architecture(param)
    , sall(nb_loc())
    , init(5)
    , so(vector<double>(0))
    , timesteps(4)
    , basal(0.1)
    , connectivity(0.1)
{
	for (int n=0; n<nb_loc(); n++)
	{
		so.push_back(init_value());
	}	
}


// functions

double ArchiRegulatory::init_value() const
{
	return init;
}


static unsigned int count_init = 0;

vector<double> ArchiRegulatory::init_pattern() const
{
	count_init++;
	vector<double> allele;
	
	for (int n=0 ; n<sall ; n++)
	{
		double value = Random::randnum();
		if (value < connectivity)
		{
			allele.push_back(1);
		}
		else
		{
			allele.push_back(0);
		}
	}
		
	std::cout << count_init << "\n";
	return(allele);
}


Phenotype ArchiRegulatory::phenotypic_value (const Genotype& genotype) const
{
	ParameterSet param;
	
	//~ vector<vector<double> > matrix;
	//~ for (int i=0 ; i<nb_loc() ; i++)
	//~ {
		//~ matrix.push_back(init_pattern());
	//~ }
	
	vector<shared_ptr<Allele> > pmatrix;
	for (int i=0 ; i<nb_loc() ; i++)
	{
		pmatrix.push_back(Allele(init_pattern()));
	}


	Haplotype haplotype;
	for (int i=0 ; i<(param.getpar(INIT_PSIZE)->GetInt()) ; i++)
	{
		for (int g=0 ; g<2 ; g++)
		{
			for (int n=0 ; n<nb_loc() ; n++)
			{
				haplotype.push_back(pmatrix[n]);
			} 
		}
	}
	
	return (haplotype[1]);
	
	
	
	
	
	
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
