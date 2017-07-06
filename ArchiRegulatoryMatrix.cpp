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

using namespace std;

BOOST_CLASS_EXPORT(ArchiRegulatoryMatrix)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(ArchiRegulatoryMatrix)

BOOST_CLASS_EXPORT(ArchiWagner)
BOOST_CLASS_EXPORT(ArchiSiegal)
BOOST_CLASS_EXPORT(ArchiM2)

std::vector<double> naive_prod(const std::vector<std::vector<double>> & m, const std::vector<double> & s) 
{
    //Expects(s.size() > 0);
    //Expects(s.size() == m.size());
    //Expects(m.size() == m[0].size());
    
    unsigned int ssz = s.size();
    
    std::vector<double> ans(ssz);
    
    for (unsigned int i = 0; i < ssz; i++) {
        ans[i] = 0.0;
        for (unsigned int j = 0; j < ssz; j++) {
            ans[i] += m[i][j] * s[j];
        }
    }
    return(ans);
}

std::vector<double> naive_prod(const std::vector<double> & s, const std::vector<std::vector<double>> & m) 
{
    //Expects(s.size() > 0);
    //Expects(s.size() == m.size());
    //Expects(m.size() == m[0].size());
    
    unsigned int ssz = s.size();    
    
    std::vector<double> ans(ssz);
    
    for (unsigned int i = 0; i < ssz; i++) {
        ans[i] = 0.0;
        for (unsigned int j = 0; j < ssz; j++) {
            ans[i] += m[j][i] * s[j];
        }
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

Phenotype ArchiRegulatoryMatrix::phenotypic_value (const Genotype& genotype, bool envir, const EpigeneticInfo & epi, bool sdinittest, bool sddynamtest) const
	// the envir parameter is useless. 
{
    // flags for optimization
    const bool run_recur = (recur != 0.0);
    const bool run_plasticity = !std::all_of(plasticity_strength.begin(), plasticity_strength.end(), [](int i) 
        { return i==0.; });
    const bool run_enviro = (Environment::dynam_disturb(sddynamtest) == 0.); // not elegant, but should do the job.
    
    
	// creation of the W matrix;
	std::vector<std::vector<double>> matrix;
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		matrix.push_back(genotype.combine_at_loc(loc, &Allele::combine_mean));
	}

    unsigned int nloc = nb_loc();
	std::vector<double> St(nloc);
	
	// Initial state of the network
	for (unsigned int i = 0 ; i<nloc ; i++)
	{
		// Warning: the following does not really make sense 
		// with binary network models!!
		// Here comes the epigenetic transmission
		if (epi.is_defined()) { 
			St[i] = epi.get_epigenet()*epi.get_phenotype()[i] + (1.-epi.get_epigenet())*So[i];
		} else {
			St[i] = So[i];
		}      		
		St[i] += Environment::init_disturb(sdinittest);
	}
    if (run_plasticity) 
        this->plasticity_v(St);
	this->haircut_v(St);

	
	// dynamic simulation
	std::vector<double> h(nloc);
	std::vector<std::vector<double>> unstability = vector<vector<double>>(calcsteps, vector<double>(nloc, 0.0));
	
	for (unsigned int t=0 ; t<timesteps ;t++)
	{
		h = naive_prod(matrix, St);
        
        this->sigma_v(h);
        
        if (run_recur)
            this->recur_v(h, St);
        
        if (run_enviro) 
            this->enviro_v(h, sddynamtest);
        
        if (run_plasticity)
            this->plasticity_v(h);
            
		this->haircut_v(h);
		
		if (t >= (timesteps-calcsteps)) 
		{
			unstability[t - (timesteps-calcsteps)] = St;
		}
	}
	
	// output (from ublas to std)
	FastIMStat stat_steps(unstability);
	vector<double> Sf_mean = stat_steps.means();
	vector<double> Sf_var = stat_steps.vars();
	
	for (unsigned int i = 0; i < Sf_mean.size(); i++)
		Sf_mean[i] += Environment::final_disturb();
	this->haircut_v(Sf_mean);
	
	return Phenotype(Sf_mean, Sf_var);
}


// Protected functions

/* Theoretically useless. Just in case, running the model on the base class just calls the identity function (no sigmoid)*/
void ArchiRegulatoryMatrix::sigma(double& h) const 
{
}

void ArchiRegulatoryMatrix::sigma_v(vector<double> & vh) const
{
    for (auto &i: vh) 
        this->sigma(i);
}

void ArchiRegulatoryMatrix::haircut(double & d) const
{
}

void ArchiRegulatoryMatrix::haircut_v(vector<double> & vh) const 
{
    for (auto &i: vh)
		this->haircut(i);
}

void ArchiRegulatoryMatrix::plasticity_v(vector<double> & vh) const
{
    for (unsigned int i = 0; i < nloc; i++) {
        vh[i] += plasticity_strength[i]*(plasticity_signal[i]-vh[i]);
    }
}

void ArchiRegulatoryMatrix::recur_v(vector<double> & vh, const std::vector<double> & oldvh) const
{
    for (unsigned int i = 0; i < nloc; i++) {
        vh[i] =  recur*oldvh[i] + (1.-recur)*vh[i];
    }
}

void ArchiRegulatoryMatrix::enviro_v(vector<double> & vh, bool sddynamtest) const
{
    for (unsigned int i = 0; i < vh.size(); i++) {
        vh[i] += Environment::dynam_disturb(sddynamtest);
    }
}


unsigned int ArchiRegulatoryMatrix::nb_phen() const
{
    return nb_loc();
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


void ArchiWagner::sigma(double & h) const 
{
	if (h < 0.) {
		h = -1.;
	} else if (h > 0.) {
		h = 1.;
	} else {
        h = 0.;
    }
}	

void ArchiWagner::sigma_v(vector<double> & vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiWagner::sigma(vh[i]);
    }
}

void ArchiWagner::haircut(double & d) const 
{
	if (d < -1.) d = -1.;
	else if (d > 1.) d = 1.;
}

void ArchiWagner::haircut_v(vector<double> & vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiWagner::haircut(vd[i]);
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


void ArchiSiegal::sigma(double &h) const 
{
	h = (2. / (1. + exp(-basal*h)) ) -1.;
}

void ArchiSiegal::sigma_v(vector<double>& vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiSiegal::sigma(vh[i]);
    }
}

void ArchiSiegal::haircut(double & d) const 
{
	if (d < -1.) d = -1.;
	else if (d > 1.) d = 1.;
}

void ArchiSiegal::haircut_v(vector<double> & vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiSiegal::haircut(vd[i]);
    }
}


ArchiM2::ArchiM2(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
	, basal(param.getpar(INIT_BASAL)->GetDouble())
{
    b1 = 1./basal - 1.;
    b2 = 1./(basal*(basal-1.))    ;
    
	double min = 0;
	double max = 1;
	
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

void ArchiM2::sigma(double& h) const 
{
	// return (1. / (1. + exp((-h/(basal*(1.-basal)))+log(1./basal-1.)) ));
    h = 1./(1.+b1*exp(b2*h));
}

void ArchiM2::sigma_v(vector<double>& vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiM2::sigma(vh[i]);
    }
}

void ArchiM2::haircut(double & d) const 
{
	if (d < 0.) d = 0.;
	else if (d > 1.) d = 1.;
}

void ArchiM2::haircut_v(vector<double>& vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiM2::haircut(vd[i]);
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
