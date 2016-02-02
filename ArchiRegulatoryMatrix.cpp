// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
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
#include "Environment.h"

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

BOOST_CLASS_EXPORT(ArchiRegulatoryMatrix)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(ArchiRegulatoryMatrix)

BOOST_CLASS_EXPORT(ArchiWagner)
BOOST_CLASS_EXPORT(ArchiSiegal)
BOOST_CLASS_EXPORT(ArchiM2)

// convert a boost vector into a regulat std::vector and vice-versa
template<class T>
std::vector<T> boost_to_std_vector(const boost::numeric::ublas::vector<T> & x) 
{
	// not very fast, optimize with std::copy if necessary
	std::vector<T> ans;
	for (unsigned int i = 0; i < x.size(); i++) 
	{
		ans.push_back(x(i));
	}
	return(ans);
}


///////////////////// Class ArchiRegulatoryMatrix //////////////////////////////////

// constructors and destuctor

/* constructor using the paramater given in Architecture and the parameters files */
ArchiRegulatoryMatrix::ArchiRegulatoryMatrix(const ParameterSet& param) 
    : Architecture(param)
    , sall(nb_loc())
    , recur(param.getpar(INIT_RECURRENCE)-> GetDouble())
    , timesteps(param.getpar(DEV_TIMESTEPS)->GetInt())
    , calcsteps(param.getpar(DEV_CALCSTEPS) -> GetInt())    
{
	init_connectivity_matrix(param); // creates connectivity_matrix
}

ArchiRegulatoryMatrix::~ArchiRegulatoryMatrix()
{
}


// functions

/* initialization of the alleles (1 allele = 1 vector of values) 
 * depends on the connectivity matrix */
shared_ptr<Allele> ArchiRegulatoryMatrix::allele_init(const ParameterSet & param, unsigned int loc) const
{
	bool clonal = (param.getpar(INIT_CLONAL)->GetString() == CL_clonal);
	vector<double> temp_allele;
	for (unsigned int i = 0; i < sall; i++) 
	{
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
	
	string type_alleles = param.getpar(TYPE_ALLELES) -> GetString();
    shared_ptr<Allele> a;
    if (type_alleles==TA_norm)
    {
        a = shared_ptr<Allele>(new Allele(temp_allele));
    }
    else if (type_alleles==TA_zero)
    {
		a = shared_ptr<Allele>(new Allele_zero(temp_allele));
    }
    else 
    {
		cerr << "Unknown allele type -- theoretically this should not happen" << endl;
		exit(EXIT_FAILURE);
	}
    return(a);
}

Phenotype ArchiRegulatoryMatrix::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi) const
{
	// creation of the W matrix;
	std::vector<std::vector<double>> matrix;
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		matrix.push_back(Allele::combine_mean(*genotype.gam_father.haplotype[loc], *genotype.gam_mother.haplotype[loc]));
	}
	
	// creation of the w_matrix and st_vector (from std to ublas)
	unsigned int nloc = nb_loc();
	boost::numeric::ublas::vector<double> St(nloc);
	boost::numeric::ublas::matrix<double> W(nloc, nloc); 
	
	// Initial state of the network
	for (unsigned int i = 0 ; i<nloc ; i++)
	{
		// Warning: the following does not really make sense 
		// with binary network models!!
		// Here comes the epigenetic transmission
		if (epi.is_defined()) { 
			St(i) = epi.get_epigenet()*epi.get_phenotype()[i] + (1.-epi.get_epigenet())*So[i];
		}
		
		if (envir) {
			St(i) += Environment::init_disturb();
		}
	}
	haircut(St);
	
	for (unsigned int i=0; i<nloc; i++) 
	{
        for (unsigned int j=0; j<nloc; j++) 
        {
            W(i,j) = matrix[i][j];            
        }
	}
	
	// dynamic simulation
	boost::numeric::ublas::vector<double> h(nloc);
	std::vector<std::vector<double>> unstability;
	
	for (unsigned int t=0 ; t<timesteps ;t++)
	{
		h = prod(St,W);
		
		for (unsigned int i=0 ; i<h.size() ; i++)
		{
			St(i) = recur*St(i) + (1-recur)*(this->sigma(h(i)));
			/* WARNING : the recurrence parameter should be put at 0 for the Wagner and Siegal model.
			 * It has no sense for the wagner and siegal model, only for the M2 model.
			 * (But implementing it only for the M2 model will be difficult due to the structure of the program) */
			//St(i) = this->sigma(h(i));
			St(i) += Environment::dynam_disturb();
		}
		haircut(St);
		
		if (t > (timesteps-calcsteps)) 
		{
			unstability.push_back(boost_to_std_vector(St));
		}
	}
	
	// output (from ublas to std)
	InvertedMStat stat_steps(unstability);
	vector<double> Sf_mean = stat_steps.means();
	vector<double> Sf_var = stat_steps.vars();
	
	for (unsigned int i = 0; i < Sf_mean.size(); i++)
		Sf_mean[i] += Environment::final_disturb();
	haircut(Sf_mean);
	
	return Phenotype(Sf_mean, Sf_var);
}


// Protected functions

/* Theoretically useless. Just in case, running the model on the base class just calls the identity function (no sigmoid)*/
double ArchiRegulatoryMatrix::sigma(double h) const 
{
	return(h);
}

//static unsigned int count_init = 0; // debug instruction

/* initialization of the connectivity matrix depend off the connectivity value */
void ArchiRegulatoryMatrix::init_connectivity_matrix(const ParameterSet & param)
{	
	double connectivity = param.getpar(INIT_CONNECT)-> GetDouble();
	double connectivity_diag = param.getpar(INIT_CONDIAG)->GetDouble();
	
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		vector<double> allele_pattern;
		for (unsigned int n=0 ; n<sall ; n++)
		{
			double threshold = connectivity;
			if (n == loc)
			{
				threshold = connectivity_diag;
			}
			
			if (Random::randnum() < threshold)
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
}

void ArchiRegulatoryMatrix::haircut(boost::numeric::ublas::vector<double> & vec) const 
{
}

void ArchiRegulatoryMatrix::haircut(vector<double> & vec) const 
{
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
		cerr << "Median So has no real significance in Wagner model." << endl << endl;
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
		cerr <<	"Basal So cannot be used in Wagner model : switch to random binary !" << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			double s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
}

ArchiWagner::~ArchiWagner()
{
	// iofile is not empty when the param constructor have been called 
	// and the FILE_ARCHI option was provided.
	if (iofile != "") {
		ofstream os(iofile);
		boost::archive::text_oarchive oa(os);
		Architecture * tmp = this;
		oa << tmp;
		// streams closed along with the destructors
	}
}


double ArchiWagner::sigma(double h) const 
{
	if (h<0)
	{
		return (-1.);
	} 
	else if (h>0) 
	{
		return (1.);
	} 
	else 
	{
		return (0.);
	}
}	

void ArchiWagner::haircut(boost::numeric::ublas::vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec(i) < -1.) vec(i) = -1.;
		else if (vec(i) > 1.) vec(i) = 1.;
	}
}

void ArchiWagner::haircut(vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec[i] < -1.) vec[i] = -1.;
		else if (vec[i] > 1.) vec[i] = 1.;
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

ArchiSiegal::~ArchiSiegal()
{
	// iofile is not empty when the param constructor have been called 
	// and the FILE_ARCHI option was provided.
	if (iofile != "") {
		ofstream os(iofile);
		boost::archive::text_oarchive oa(os);
		Architecture * tmp = this;
		oa << tmp;
		// streams closed along with the destructors
	}
}


double ArchiSiegal::sigma(double h) const 
{
	return ((2. / (1. + exp(-basal*h)) ) -1.);
}

void ArchiSiegal::haircut(boost::numeric::ublas::vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec(i) < -1.) vec(i) = -1.;
		else if (vec(i) > 1.) vec(i) = 1.;
	}
}

void ArchiSiegal::haircut(std::vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec[i] < -1.) vec[i] = -1.;
		else if (vec[i] > 1.) vec[i] = 1.;
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

double ArchiM2::sigma(double h) const 
{
	return (1. / (1. + exp((-h/(basal*(1.-basal)))+log(1./basal-1.)) ));
}

void ArchiM2::haircut(boost::numeric::ublas::vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec(i) < 0.) vec(i) = 0.;
		else if (vec(i) > 1.) vec(i) = 1.;
	}
}

void ArchiM2::haircut(std::vector<double> & vec) const 
{
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec[i] < 0.) vec[i] = 0.;
		else if (vec[i] > 1.) vec[i] = 1.;
	}
}

ArchiM2::~ArchiM2()
{
	// iofile is not empty when the param constructor have been called 
	// and the FILE_ARCHI option was provided.
	if (iofile != "") {
		ofstream os(iofile);
		boost::archive::text_oarchive oa(os);
		Architecture * tmp = this;
		oa << tmp;
		// streams closed along with the destructors
	}
}
