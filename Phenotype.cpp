// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Phenotype.h"

#include <cassert>

using namespace std;



// constructors and destructors

Phenotype::Phenotype()
{
	// Nothing yet.
}


Phenotype::Phenotype(double init) 
{
	// only one phenotypic dimension!
	pheno.push_back(init);
}


Phenotype::Phenotype(const vector<double> & init) 
{
	pheno = init;
}


Phenotype::Phenotype(const Phenotype & templ)
{
	pheno = templ.pheno;
}


Phenotype::~Phenotype() 
{
	
}


// operator overload 

Phenotype & Phenotype::operator = (const Phenotype & templ) 
{
	if (this == &templ)
        return (*this);
    pheno = templ.pheno;
    return(*this);
}


// functions

double Phenotype::operator[] (const unsigned int index) const 
{
	assert(index < dimensionality());
	return(pheno[index]);
}


unsigned int Phenotype::dimensionality() const 
{
	return(pheno.size());
}


void Phenotype::write_debug (ostream& out) const
{
	for (unsigned int i = 0; i < dimensionality(); i++) {
		out << pheno[i];
		if (i < dimensionality()-2)
			out << "/";
	}
	out << endl;
}


void Phenotype::write_simple (ostream& out) const
{
	write_debug(out);
}

ostream& operator << (ostream& out, const Phenotype& phen)
{
	for (unsigned int k = 0; k < phen.dimensionality(); k++) {
		out << phen[k] << "\t";
	}
	return(out);
}

/////////////////////// class PhenotypeStat

PhenotypeStat::PhenotypeStat(const vector<Phenotype> & vec_phen)
	: MultivariateStat(transpose_phen_matrix(vec_phen))
{
	
}

// static function called into the initialization list of the constructor
vector<vector<double> > PhenotypeStat::transpose_phen_matrix(const vector<Phenotype> & vec_phen)
{ // the vector of vectors has to be transposed in order to get means and variances of traits
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
