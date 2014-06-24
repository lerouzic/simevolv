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



#include "ArchiRegulatoryMatrix.h"
#include "Parconst.h"
#include "Random.h"
#include "Statistics.h"
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


// convert a boost vector into a regulat std::vector and vice-versa

template<class T>
std::vector<T> boost_to_std(const boost::numeric::ublas::vector<T> & x) 
{
	// not very fast, optimize with std::copy if necessary
	std::vector<T> ans;
	for (unsigned int i = 0; i < x.size(); i++) 
	{
		ans.push_back(x(i));
	}
	return(ans);
}


// constructors and destuctor

/* default constructor  -  should never be used */
ArchiRegulatoryMatrix::ArchiRegulatoryMatrix()
{
    assert(false); // The default constructor should never be used.
}


/* constructor using the paramater given in Architecture and the parameters files */
ArchiRegulatoryMatrix::ArchiRegulatoryMatrix(const ParameterSet& param) 
    : Architecture(param)
    , sall(nb_loc())
    , timesteps(param.getpar(DEV_TIMESTEPS)->GetInt())
    , calcsteps(param.getpar(DEV_CALCSTEPS) -> GetInt())
{
	init_connectivity_matrix(param); // creates connectivity_matrix
}


// functions

/* initialization of the alleles (1 allele = 1 vector of value) 
 * depend on the connectivity matrix */
shared_ptr<Allele> ArchiRegulatoryMatrix::allele_init(const ParameterSet & param, unsigned int loc) const
{
	bool clonal = (param.getpar(INIT_CLONAL)->GetString() == CL_clonal);
	vector<double> temp_allele;
	for (unsigned int i = 0; i < sall; i++) {
		assert(connectivity_matrix.size() > loc);
		assert(connectivity_matrix[loc].size() > i);
		if (connectivity_matrix[loc][i] == 0.0) 
		{
			temp_allele.push_back(0.0);
		} 
		else 
		{
			double value;
			if (clonal) 
			{
				value = connectivity_matrix[loc][i];
			} 
			else 
			{
				value = param.getpar(INIT_ALLELES) -> GetDouble();
			}
			temp_allele.push_back(value);
		}
	}
	shared_ptr<Allele> a (new Allele(temp_allele)); // Here the class of Allele should be determined. 
	return(a);
}


//static unsigned int count_init = 0; // debug instruction

/* initialization of the connectivity matrix
 * depend off tne connectivity value */
void ArchiRegulatoryMatrix::init_connectivity_matrix(const ParameterSet & param)
{
	//count_init++; // debug instruction;
	
	double connectivity = param.getpar(INIT_CONNECT)-> GetDouble();
	
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		vector<double> allele_pattern;
		for (unsigned int n=0 ; n<sall ; n++)
		{
			if (Random::randnum() < connectivity)
			{
				/* This has some interest only in clonal populations. In non-clonal pops, this value will be overwritten anyway. 
				 * double value = 1.0  - in non-clonal, this would be enough */
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
	//std::cout << count_init << "\n";  // debug instruction
}


Phenotype ArchiRegulatoryMatrix::phenotypic_value (const Genotype& genotype) const
{
	// creation of the W matrix;
	std::vector<std::vector<double> > matrix;
	for (unsigned int loc=0; loc < nb_loc(); loc++) 
	{
		matrix.push_back(Allele::combine_add(*genotype.gam_father.haplotype[loc], *genotype.gam_mother.haplotype[loc]));
	}
	
	// creation of the w_matrix and st_vector (from std to ublas)
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
	std::vector<std::vector<double> > unstability;
	
	for (unsigned int t=0 ; t<timesteps ;t++)
	{
		h = prod(St,W);
		for (unsigned int i=0 ; i<h.size() ; i++)
		{
			St(i) = sigma(h(i));
			//~ cout << "Phen " << i << " time " << t << ": " << St(i) << "\n";
		}
		if (t > (timesteps-calcsteps)) 
		{
			unstability.push_back(boost_to_std(St));
		}
	}

	// output (from ublas to std)
	InvertedMStat stat_steps(unstability);
	vector<double> Sf_mean = stat_steps.means();
	vector<double> Sf_var = stat_steps.vars();
	
	return Phenotype(Sf_mean, Sf_var);
}




//////////////////////////// INHERITED CLASSES //////////////////////////////////

// constructors 

ArchiWagner::ArchiWagner(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
{ 
	int min = -1;
	int max = 1;
	
	string type_so = param.getpar(TYPE_SO)->GetString();
    if (type_so==SO_min)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(min);
		}	
    }
    else if (type_so==SO_max)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(max);
		}	
	}
	else if (type_so==SO_med)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);			
		}	
	}
	else if (type_so==SO_randbin)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		cerr << "Random So has no real significance in Wagner model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = (max-min)*Random::randnum()+min;
			So.push_back(s);
		}
	}
	else
	{
		cerr <<	"Basal So cannot be used in Wagner model : switch to Median So." << endl;
		cerr << "Median So has no real significance in Wagner model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);
		}	
	}
}


ArchiMasel::ArchiMasel(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
{ 
	int min = 0;
	int max = 1;
	
	string type_so = param.getpar(TYPE_SO)->GetString();
    if (type_so==SO_min)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(min);
		}	
    }
    else if (type_so==SO_max)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(max);
		}	
	}
	else if (type_so==SO_med)
	{
		cerr << "Median So has no real significance in Masel model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);	
		}	
	}
	else if (type_so==SO_randbin)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		cerr << "Random So has no real significance in Masel model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = (max-min)*Random::randnum()+min;
			So.push_back(s);
		}
	}
	else 
	{
		cerr <<	"Basal So cannot be used in Masel model : switch to Median So." << endl;
		cerr << "Median So has no real significance in Masel model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);
		}	
	}	
}


ArchiSiegal::ArchiSiegal(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
	, basal(param.getpar(INIT_BASAL)->GetDouble())
{ 
	int min = -1;
	int max = 1;
	
	string type_so = param.getpar(TYPE_SO)->GetString();
    if (type_so==SO_min)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(min);
		}	
    }
    else if (type_so==SO_max)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(max);
		}	
	}
	else if (type_so==SO_med)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);			
		}	
	}
	else if (type_so==SO_randbin)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = (max-min)*Random::randnum()+min;
			So.push_back(s);
		}
	}
	else
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(basal);
		}
	}
}


ArchiM2::ArchiM2(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
	, basal(param.getpar(INIT_BASAL)->GetDouble())
{ 
	int min = 0;
	int max = 1;
	
	string type_so = param.getpar(TYPE_SO)->GetString();
    if (type_so==SO_min)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(min);
		}	
    }
    else if (type_so==SO_max)
    {
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(max);
		}	
	}
	else if (type_so==SO_med)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back((min+max)/2);			
		}	
	}
	else if (type_so==SO_randbin)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = (max-min)*Random::randnum()+min;
			So.push_back(s);
		}
	}
	else
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			So.push_back(basal);
		}
	}
}
