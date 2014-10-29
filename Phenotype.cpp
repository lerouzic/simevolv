// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
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



#include "Phenotype.h"
#include "OutputFormat.h"

#include <cassert>

using namespace std;

ostream& operator << (ostream& out, const Phenovec & pvec)
{
	for (unsigned int i = 0; i < pvec.size(); i++) {
		out << pvec[i];
		if (i < pvec.size() - 1)
			out << "\t";
	}
	return(out);
}

void outformat(std::ostream & out, const Phenovec & pv, unsigned int width /*=10*/, 
			unsigned int precision /* =5 */, std::string sep /* ="" */) {
	for (unsigned int i = 0; i < pv.size(); i++)
		outformat(out, pv[i], width, precision, sep);
}

// constructors and destructors

/* default constructor */
Phenotype::Phenotype()
{
}


/* constructor using an initialization value 
 * for uni-dimension phentotype */
Phenotype::Phenotype(double init) 
{
	pheno.push_back(init);
	unstabpheno.push_back(0.0);
}


/* constructor using a vector of initialization value 
 * for multi-dimension phenotype */
Phenotype::Phenotype(const vector<double> & init) 
{
	pheno = init;
	// when the unstability vector is not provided, we assume that phenotypes
	// are all stable.
	for (unsigned int i = 0; i < pheno.size(); i++) 
	{
		unstabpheno.push_back(0.0);
	}
}

Phenotype::Phenotype(const Phenovec & init)
{
	pheno = init;
	// when the unstability vector is not provided, we assume that phenotypes
	// are all stable.
	for (unsigned int i = 0; i < pheno.size(); i++) 
	{
		unstabpheno.push_back(0.0);
	}
}

Phenotype::Phenotype(const vector<double>& init, const vector<double>& initunstab)
{
	pheno = init;
	unstabpheno = initunstab;
}


/* constructor using a previous phenotype as pattern */
Phenotype::Phenotype(const Phenotype & templ)
{
	pheno = templ.pheno;
	unstabpheno = templ.unstabpheno;
}


/* destructor */
Phenotype::~Phenotype() 
{
	
}


// operator overload 

Phenotype & Phenotype::operator = (const Phenotype & templ) 
{
	if (this == &templ)
        return (*this);
    pheno = templ.pheno;
    unstabpheno = templ.unstabpheno;
    return(*this);
}

// Warning: this returns the phenotypic value!
double Phenotype::operator[] (const unsigned int index) const 
{
	return(get_pheno(index));
}

Phenovec Phenotype::get_pheno() const
{
	return(pheno);
}

double Phenotype::get_pheno(const unsigned int index) const
{
	assert(index < dimensionality());
	return(pheno[index]);
}

Phenovec Phenotype::get_unstab() const
{
	return(unstabpheno);
}

double Phenotype::get_unstab(const unsigned int index) const
{
	assert(index < dimensionality());
	assert(pheno.size() == unstabpheno.size());
	return(unstabpheno[index]);
}

void Phenotype::add_pheno(const unsigned int index, const double effect) 
{
	assert (index < dimensionality());
	pheno[index] += effect;
}


// functions

/* return the dimensionnality of the phenotype (number of phenotypes observed) */
unsigned int Phenotype::dimensionality() const 
{
	assert(pheno.size() == unstabpheno.size());
	return(pheno.size());
}


void Phenotype::write_debug (ostream& out) const
{
	for (unsigned int i = 0; i < dimensionality(); i++) 
	{
		//out << pheno[i] << "(" << unstabpheno[i] << ")"; 	//for debug
		out << pheno[i];
		if (i < dimensionality()-2)
			out << "/";
	}
	out << endl;
}


// Output

void Phenotype::write_simple (ostream& out) const
{
	write_debug(out);
}

ostream& operator << (ostream& out, const Phenotype& phen)
{
	for (unsigned int k = 0; k < phen.dimensionality(); k++) 
	{
		//out << phen.pheno[k] << "(" << phen.unstabpheno[k] << ")\t";   // for debug
		out << phen.pheno[k] << "\t";
	}
	return(out);
}



void outformat(std::ostream & out, const Phenotype & ph, unsigned int width /*=10*/, 
	unsigned int precision /*=5*/, std::string sep /*=""*/) {
		outformat(out, ph.pheno, width, precision, sep);
}



/////////////////////// class PhenotypeStat

// constructors and destructors

/* constructor using a vector of phenotype */
PhenotypeStat::PhenotypeStat(const vector<Phenotype> & vec_phen)
{
	assert(vec_phen.size() > 2);
	pheno = new MultivariateStat(transpose_phen_matrix(vec_phen));
	unstab = new MultivariateStat(transpose_unstabphen_matrix(vec_phen));
}

PhenotypeStat::PhenotypeStat(const vector<Phenovec> & vec_phen)
{
	assert(vec_phen.size() > 2);
	pheno = new MultivariateStat(transpose_phenovec_matrix(vec_phen));
	unstab = NULL;
}

PhenotypeStat::~PhenotypeStat()
{
	// Very important; these are raw pointers and memory MUST be cleared when the object is destroyed!
	delete pheno;
	delete unstab; 
}

// functions

/* transpose the phen_matrix (vector of vector) to get means and variances of traits
 * static function called into the initialization list of the constructor */
vector<vector<double> > PhenotypeStat::transpose_phen_matrix(const vector<Phenotype> & vec_phen)
{ 
	assert(!vec_phen.empty());
	
	unsigned int dim = vec_phen[0].dimensionality();
	vector<vector<double> > ans;
	
	for (unsigned int k = 0; k < dim; k++) {
		vector<double> tmp;
		for (unsigned int i = 0; i < vec_phen.size(); i++) {
			tmp.push_back(vec_phen[i].get_pheno(k));
		}
		ans.push_back(tmp);
	}
	return(ans);
}

vector<vector<double> > PhenotypeStat::transpose_unstabphen_matrix(const vector<Phenotype> & vec_phen)
{ 
	assert(!vec_phen.empty());
	
	unsigned int dim = vec_phen[0].dimensionality();
	vector<vector<double> > ans;
	
	for (unsigned int k = 0; k < dim; k++) {
		vector<double> tmp;
		for (unsigned int i = 0; i < vec_phen.size(); i++) {
			tmp.push_back(vec_phen[i].get_unstab(k));
		}
		ans.push_back(tmp);
	}
	return(ans);
}

vector<vector<double> > PhenotypeStat::transpose_phenovec_matrix(const vector<Phenovec> & vec_phen)
{
	assert(!vec_phen.empty());
	
	unsigned int dim = vec_phen[0].dimensionality();
	vector<vector<double> > ans;
	
	for (unsigned int k = 0; k < dim; k++) {
		vector<double> tmp;
		for (unsigned int i = 0; i < vec_phen.size(); i++) {
			tmp.push_back(vec_phen[i][k]);
		}
		ans.push_back(tmp);
	}
	return(ans);
}

// Inverted multivariate stats (to be copied in multivariateStats?)

InvertedMStat::InvertedMStat(const vector<vector<double> > & vec_vec_d)
	: MultivariateStat(transpose_double_matrix(vec_vec_d))
{
}


// Exactly the same code as above (:-@!!! )
vector<vector<double> > InvertedMStat::transpose_double_matrix(const vector<vector<double> > & vec_vec_d) 
{
	assert(!vec_vec_d.empty());
	
	unsigned int dim = vec_vec_d[0].size();
	vector<vector<double> > ans;
	for (unsigned int k = 0; k < dim; k++) {
		vector<double> tmp;
		for (unsigned int i = 0; i < vec_vec_d.size(); i++) {
			if (k==0) assert(vec_vec_d[i].size() == dim);
			tmp.push_back(vec_vec_d[i][k]);
		}
		ans.push_back(tmp);
	}
	return(ans);
}
