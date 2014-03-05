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



#include "ArchiRegulatoryWagner.h"
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

ArchiRegulatoryWagner::ArchiRegulatoryWagner()
{
    assert(false); // The default constructor should never be used.
}


ArchiRegulatoryWagner::ArchiRegulatoryWagner(const ParameterSet& param) 
    : Architecture(param)
    , sall(nb_loc())
    , timesteps(param.getpar(DEV_TIMESTEPS)->GetDouble())
    , basal(param.getpar(INIT_BASAL)->GetDouble())
{
	init_connectivity_matrix(param); // creates connectivity_matrix

	for (unsigned int n=0; n<nb_loc(); n++)
	{
		So.push_back(basal);
	}	
}


// functions

shared_ptr<Allele> ArchiRegulatoryWagner::allele_init(const ParameterSet & param, unsigned int loc) const
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


//static unsigned int count_init = 0; // debug instruction


void ArchiRegulatoryWagner::init_connectivity_matrix(const ParameterSet & param)
{
	//count_init++; // debug instruction;
	
	double connectivity = param.getpar(INIT_CONNECT)->GetDouble();
	//bool clonal = (param.getpar(INIT_CLONAL)->GetString() == CL_clonal);
	
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		vector<double> allele_pattern;
		for (unsigned int n=0 ; n<sall ; n++)
		{
			if (Random::randnum() < connectivity)
			{
				// This has some interest only in clonal populations. In non-clonal pops, this value will be overwritten anyway. 
				// double value = 1.0 		// in non-clonal, this would be enough.
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
	//std::cout << count_init << "\n";
}

//~ 
//~ Phenotype ArchiRegulatoryWagner::phenotypic_value (const Genotype& genotype) const
//~ {
	//~ // creation of the W matrix;
	//~ vector<vector<double> > matrix;
	//~ for (unsigned int loc=0; loc < nb_loc(); loc++) 
	//~ {
		//~ matrix.push_back(Allele::combine_add(*genotype.gam_father.haplotype[loc], *genotype.gam_mother.haplotype[loc]));
	//~ }
	//~ 
	//~ // creation of the w_matrix and st_vector for using ublas 
	//~ //using namespace boost::numeric::ublas;
	//~ 
	//~ unsigned int nloc = nb_loc();
	//~ boost::numeric::ublas::vector<double> St(nloc);
	//~ boost::numeric::ublas::matrix<double> W(nloc, nloc); 
		//~ 
	//~ for (unsigned int i=0 ; i<nloc ; i++)
	//~ {
		//~ St(i) = So[i];
	//~ }
	//~ 
	//~ for (unsigned int i=0; i<nloc; i++) 
	//~ {
        //~ for (unsigned int j=0; j<nloc; j++) 
        //~ {
            //~ W(i,j) = matrix[i][j];
        //~ }
	//~ }
		//~ 
	//~ // simulation
	//~ boost::numeric::ublas::vector<double> h(nloc);
	//~ 
	//~ for (unsigned int t=0 ; t<timesteps ;t++)
	//~ {
		//~ h = prod(St,W);
		//~ for (unsigned int i=0 ; i<h.size() ; i++)
		//~ {
			//~ if (h(i)<0) 
			//~ {
				//~ St(i) = -1.;
			//~ }
			//~ else if (h(i)>0) 
			//~ {
				//~ St(i) = 1.;
			//~ }
			//~ else 
			//~ {
				//~ St(i) = 0.;
			//~ }
		//~ }
	//~ }
	//~ 
	//~ // output
	//~ std::vector<double> Sf;
	//~ for (unsigned int i=0 ; i<nloc ; i++)
	//~ {
		//~ Sf.push_back(St(i));
	//~ }
	//~ 
	//~ return Phenotype(Sf);
//~ }
