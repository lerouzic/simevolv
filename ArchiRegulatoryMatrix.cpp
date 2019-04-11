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

#ifdef SERIALIZATION_TEXT
BOOST_CLASS_EXPORT(ArchiRegulatoryMatrix)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(ArchiRegulatoryMatrix)

BOOST_CLASS_EXPORT(ArchiWagner)
BOOST_CLASS_EXPORT(ArchiSiegal)
BOOST_CLASS_EXPORT(ArchiM2)
#endif

template<typename T>
std::vector<T> naive_prod(const std::vector<std::vector<T>> & m, const std::vector<T> & s) 
{
    //Expects(s.size() > 0);
    //Expects(s.size() == m.size());
    //Expects(m.size() == m[0].size());
    
    unsigned int ssz = s.size();
    
    std::vector<T> ans(ssz);
    
    for (unsigned int i = 0; i < ssz; i++) {
        ans[i] = 0.0;
        for (unsigned int j = 0; j < ssz; j++) {
            ans[i] += m[i][j] * s[j];
        }
    }
    return(ans);
}

template<typename T>
std::vector<T> naive_prod(const std::vector<T> & s, const std::vector<std::vector<T>> & m) 
{
    //Expects(s.size() > 0);
    //Expects(s.size() == m.size());
    //Expects(m.size() == m[0].size());
    
    unsigned int ssz = s.size();    
    
    std::vector<T> ans(ssz);
    
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
    , recur(param.getpar(INIT_RECURRENCE)-> GetDouble())
    , timesteps(param.getpar(DEV_TIMESTEPS)->GetInt())
    , calcsteps(param.getpar(DEV_CALCSTEPS) -> GetInt())    
{
	// update_param_internal(param); // This should be called in derived classes only
	sall = nb_loc();
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
	vector<allele_type> temp_allele;
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
			allele_type value;
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
    const bool run_plasticity = !std::all_of(plasticity_strength.begin(), plasticity_strength.end(), [](double i) 
        { return i==0.; });
    const bool run_enviro = (Environment::dynam_disturb(sddynamtest) == 0.); // not elegant, but should do the job.
    
    
	// creation of the W matrix;
	std::vector<std::vector<allele_type>> matrix;
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		matrix.push_back(genotype.combine_at_loc(loc, &Allele::combine_mean));
	}

    unsigned int nloc = nb_loc();
	std::vector<pheno_type> St(nloc);
	
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
	std::vector<pheno_type> St_next(nloc);
	std::vector<std::vector<pheno_type>> unstability = vector<vector<pheno_type>>(calcsteps, vector<pheno_type>(nloc, 0.0));
	
	for (unsigned int t=0 ; t<timesteps ;t++)
	{
		St_next = naive_prod(matrix, St);
        
        this->sigma_v(St_next);
        
        if (run_recur)
            this->recur_v(St_next, St);
        
        if (run_enviro) 
            this->enviro_v(St_next, sddynamtest);
        
        if (run_plasticity)
            this->plasticity_v(St_next);
            
		this->haircut_v(St_next);
		
		St = St_next;
		
		if (t >= (timesteps-calcsteps)) 
		{
			unstability[t - (timesteps-calcsteps)] = St;
		}
	}
	
	FastIMStat<pheno_type> stat_steps(unstability);
	vector<pheno_type> Sf_mean = stat_steps.means();
	vector<pheno_type> Sf_var = stat_steps.vars();
	
	for (unsigned int i = 0; i < Sf_mean.size(); i++)
		Sf_mean[i] += Environment::final_disturb();
	this->haircut_v(Sf_mean);
	
	Phenotype ans(PhenoTranscriptome(Sf_mean, Sf_var));
	ans.scale_transform(transfo);
	return ans;
}


// Protected functions

/* Theoretically useless. Just in case, running the model on the base class just calls the identity function (no sigmoid)*/
void ArchiRegulatoryMatrix::sigma(pheno_type& h) const 
{
}

void ArchiRegulatoryMatrix::sigma_v(vector<pheno_type> & vh) const
{
    for (auto &i: vh) 
        this->sigma(i);
}

void ArchiRegulatoryMatrix::haircut(pheno_type & d) const
{
}

void ArchiRegulatoryMatrix::haircut_v(vector<pheno_type> & vh) const 
{
    for (auto &i: vh)
		this->haircut(i);
}

void ArchiRegulatoryMatrix::plasticity_v(vector<pheno_type> & vh) const
{
    // a bit optimized...
    auto vh_it = vh.begin();
    auto pst_it = plasticity_strength.begin(); // this is probably a const_iterator because of the const method
    auto psi_it = plasticity_signal.begin();
    for ( ; vh_it != vh.end(); vh_it++, pst_it++, psi_it++) {
        *vh_it += (*pst_it) * (*psi_it - *vh_it);
    }
}

void ArchiRegulatoryMatrix::recur_v(vector<pheno_type> & vh, const std::vector<pheno_type> & oldvh) const
{
    auto vh_it = vh.begin();
    auto oldvh_it = oldvh.begin();
    for ( ; vh_it != vh.end(); vh_it++, oldvh_it++) {
        *vh_it += recur*(*oldvh_it - *vh_it);
    }
}

void ArchiRegulatoryMatrix::enviro_v(vector<pheno_type> & vh, bool sddynamtest) const
{
    for (auto &i: vh)
        i += Environment::dynam_disturb(sddynamtest);
}


unsigned int ArchiRegulatoryMatrix::nb_phen() const
{
    return nb_loc();
}

//static unsigned int count_init = 0; // debug instruction

/* initialization of the connectivity matrix depend off the connectivity value */
void ArchiRegulatoryMatrix::init_connectivity_matrix(const ParameterSet & param)
{	
	rate_type connectivity = param.getpar(INIT_CONNECT)-> GetDouble();
	rate_type connectivity_diag = param.getpar(INIT_CONDIAG)->GetDouble();
	
	for (unsigned int loc = 0; loc < nb_loc(); loc++) 
	{
		vector<allele_type> allele_pattern;
		for (unsigned int n=0 ; n<sall ; n++)
		{
			rate_type threshold = connectivity;
			if (n == loc)
			{
				threshold = connectivity_diag;
			}
			
			if (Random::randnum() < threshold)
			{
				/* This has some interest only in clonal populations. In non-clonal pops, this value will be overwritten anyway. 
				 * double value = 1.0  - in non-clonal, this would be enough */
				allele_type value = param.getpar(INIT_ALLELES) -> GetDouble();
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
    update_param_internal(param); 
 
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
			pheno_type s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		cerr << "Random So has no real significance in Wagner model." << endl << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			pheno_type s = (max-min)*Random::randnum()+min;
			So.push_back(s);
		}
	}
	else
	{
		cerr <<	"Basal So cannot be used in Wagner model : switch to random binary !" << endl;
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			pheno_type s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
}

ArchiWagner::~ArchiWagner()
{
}

void ArchiWagner::sigma(pheno_type & h) const 
{
	if (h < 0.) {
		h = -1.;
	} else if (h > 0.) {
		h = 1.;
	} else {
        h = 0.;
    }
}	

void ArchiWagner::sigma_v(vector<pheno_type> & vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiWagner::sigma(vh[i]);
    }
}

void ArchiWagner::haircut(pheno_type & d) const 
{
	if (d < -1.) d = -1.;
	else if (d > 1.) d = 1.;
}

void ArchiWagner::haircut_v(vector<pheno_type> & vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiWagner::haircut(vd[i]);
    }
}

ArchiSiegal::ArchiSiegal(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
	, basal(param.getpar(INIT_BASAL)->GetDouble())
{
	update_param_internal(param); 
    
     
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
			pheno_type s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			pheno_type s = (max-min)*Random::randnum()+min;
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
}


void ArchiSiegal::sigma(pheno_type &h) const 
{
	h = (2. / (1. + exp(-basal*h)) ) -1.;
}

void ArchiSiegal::sigma_v(vector<pheno_type>& vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiSiegal::sigma(vh[i]);
    }
}

void ArchiSiegal::haircut(pheno_type & d) const 
{
	if (d < -1.) d = -1.;
	else if (d > 1.) d = 1.;
}

void ArchiSiegal::haircut_v(vector<pheno_type> & vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiSiegal::haircut(vd[i]);
    }
}


ArchiM2::ArchiM2(const ParameterSet& param) 
	: ArchiRegulatoryMatrix(param)
	, basal(param.getpar(INIT_BASAL)->GetDouble())
{
	update_param_internal(param); 
    
    b1 = 1./basal - 1.;
    b2 = 1./(basal*(basal-1.))    ;
    
	pheno_type min = 0;
	pheno_type max = 1;
	
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
			pheno_type s = Random::randnum();
			if (s < 0.5) { So.push_back(min); }
			else { So.push_back(max); }
		}	
	}
	else if (type_so==SO_rand)
	{
		for (unsigned int n=0; n<nb_loc(); n++)
		{
			pheno_type s = (max-min)*Random::randnum()+min;
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

void ArchiM2::sigma(pheno_type& h) const 
{
	// return (1. / (1. + exp((-h/(basal*(1.-basal)))+log(1./basal-1.)) ));
    h = 1./(1.+b1*exp(b2*h));
}

void ArchiM2::sigma_v(vector<pheno_type>& vh) const
{
    for (unsigned int i = 0; i < vh.size(); ++i) {
        ArchiM2::sigma(vh[i]);
    }
}

void ArchiM2::haircut(pheno_type& d) const 
{
	if (d < 0.) d = 0.;
	else if (d > 1.) d = 1.;
}

void ArchiM2::haircut_v(vector<pheno_type>& vd) const
{
    for (unsigned int i = 0; i < vd.size(); ++i) {
        ArchiM2::haircut(vd[i]);
    }
}

ArchiM2::~ArchiM2()
{
}
