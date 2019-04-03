// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "OutputFormat.h"

#include "Phenotype.h"

#include <sstream>

using namespace std;


class Phenotype;


void outformat(ostream & out, const double number, unsigned int width /*=12*/, unsigned int precision /*=5*/, const string & sep /*=""*/) 
{
	out << setw(width) << setprecision(precision) << left << number << sep;
}
void outformat(ostream & out, const int number, unsigned int width /*=12*/, const string & sep /*=""*/) 
{
	out << setw(width) << left << number << sep;
}

void outformat(ostream & out, const vector<double> & numbers, unsigned int width /*=12*/, unsigned int precision /*=5*/, const string & sep /*=""*/) 
{
	for (unsigned int i = 0; i < numbers.size(); i++)
		outformat(out, numbers[i], width, precision, sep);
}

void outformat(ostream & out, const string & text, unsigned int width /*=12*/, const string & sep /*=""*/) 
{
	out << setw(width) << left << text << sep;
}

void outformat(ostream & out, unsigned int index, const string & text, unsigned int width /*=12*/, const string & sep /*=""*/) 
{
	ostringstream o;
	o << text << index;
	out << setw(width) << left << o.str() << sep;
}

void outformat(ostream & out, unsigned int index_i, unsigned int index_j, const string & text, unsigned int width /*=12*/, const string & sep /*=""*/) 
{
	ostringstream o;
	o << text << index_i << "_" << index_j;
	out << setw(width) << left << o.str() << sep;
}


void outformat(ostream & out, const Phenotype & pheno,
                      unsigned int width, unsigned int precision, string sep)
{
    assert(pheno.is_defined());
    pheno.pheno_ptr->outformat(out, width, precision, sep);
}

void outformat2(ostream & out, const Phenotype & pheno,
                      unsigned int width, unsigned int precision, string sep)
{
    assert(pheno.is_defined());
    pheno.pheno_ptr->outformat2(out, width, precision, sep);
}
