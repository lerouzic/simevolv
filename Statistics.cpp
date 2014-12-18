// Copyright 2014       Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Statistics.h"

#include <cassert>
#include <iostream>
#include <cmath>

using namespace std;



////////////////////////////////// UNIVARIATE ////////////////////////////////////////

// Constructors and initialization 

UnivariateStat::UnivariateStat(const vector<double> & vv)
	: data(vv)
	, sum_i(0.0)
	, sum_i2(0.0)
{
	if (data.size() < 2) 
	{
		cerr << "A data set of size " << data.size() << " is not eligible for statistics." << endl;
		assert("Terminate.");
	}
	initialize();
}

/* computes and stores intermediate results. */
void UnivariateStat::initialize()
{
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		if (isfinite(data[i])) 
		{
			sum_i += data[i];
			sum_i2 += data[i]*data[i];
		} 
		else 
		{
			// cerr << "Warning: trying to compute statistics with non-finite numbers" << endl;
			data.erase(data.begin()+i);
		}
	}
}


// functions

double UnivariateStat::mean() const {
	return(sum_i / static_cast<double>(data.size()));
}

double UnivariateStat::var() const {
	double mm = mean();
	return(sum_i2/static_cast<double>(data.size()) - mm*mm);
}


////////////////////////////////// MULTIVARIATE /////////////////////////////////

MultivariateStat::MultivariateStat(const vector<vector<double> > & v)
	: data(v)
{ 
	// It is a bit complicated to assert here whether the vector<vector<double> > is OK. 
	// Temporary storage variables as well as checks on the dimensions will occur afterwards, 
	// in initialize(), which computes the sums of x_i and x_i^2. 
	initialize();
}

void MultivariateStat::initialize()
{
	unsigned int size1;
	
	// Tests if the data behaves properly (same number of measurements for each category)
	assert(data.size() > 0);
	size1 = data[0].size();
	assert(size1 > 0);
	if (data.size() > 1) 
	{
		for (unsigned int i = 1; i < data.size(); i++) 
		{ // i = 1: desired behaviour
			assert(data[i].size() == size1);
		}
	}
	
	// creates the matrix of sum_ij
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		vector<double> tmp_sumij(data.size(), 0.0);
		sum_ij.push_back(tmp_sumij);
	}
	
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		
		// vector of sum_i
		double tmp_sum = 0.0;
		for (unsigned int k = 0; k < data[i].size(); k++) 
		{
			tmp_sum += data[i][k];
		}
		sum_i.push_back(tmp_sum);
		
		// matrix of sum_ij
		for (unsigned int j = i; j < data.size(); j++) 
		{
			double tmp_sumij = 0.0;
			for (unsigned int k = 0; k < data[i].size(); k++)
			{
				tmp_sumij += data[i][k] * data[j][k];
			}
			sum_ij[i][j] = tmp_sumij;
			if (i != j) { // does not harm, but useless when i == j
				sum_ij[j][i] = tmp_sumij; 
			}
		}
	}
}

vector<double> MultivariateStat::means() const
{
	vector<double> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = sum_i[i]/static_cast<double>(data[i].size());
	}
	return(ans);
}

vector<double> MultivariateStat::vars() const
{
	vector<double> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = var(i);
	}
	return(ans);
}

vector<vector<double> > MultivariateStat::vcov() const
{
	vector<vector<double> > ans (data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		vector<double> tmp (data.size());
		for (unsigned int j = 0; j < data.size(); j++) 
		{
			tmp[j] = cov(i, j);
		}
		ans.push_back(tmp);
	}
	return(ans);
}

double MultivariateStat::mean(unsigned int i) const
{
	assert (i < data.size());
	return(sum_i[i] / static_cast<double>(data[i].size()));
}

double MultivariateStat::var(unsigned int i) const 
{
	assert (i < data.size());
	return(sum_ij[i][i]/static_cast<double>(data[i].size()) - mean(i)*mean(i));
}

double MultivariateStat::cov(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(sum_ij[i][j]/static_cast<double>(data[i].size()) - mean(i)*mean(j));
}

double MultivariateStat::cor(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cov(i, j)/(sqrt(var(i)*var(j))));
}

double MultivariateStat::r2(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cor(i, j)*cor(i, j));	
}

double MultivariateStat::regression_slope(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cov(i, j)/var(j));
}

ostream & operator << (ostream & os, const MultivariateStat & obj) 
{
	os << "\t";
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << "C" << cat+1 << "\t";
	}
	os << "\n";
	for (unsigned int rep = 0; rep < obj.data[0].size(); rep++) 
	{
		os << rep+1 << "\t";
		for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
		{
			os << obj.data[cat][rep] << "\t";
		}
		os << "\n";
	}
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << "-------" << "\t";
	}	
	os << "\n";
	
	os << "Sum i" << "\t";
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << obj.sum_i[cat] << "\t";
	}	
	os << "\n";
	os << "Sum i2" << "\t";
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << obj.sum_ij[cat][cat] << "\t";
	}	
	os << "\n";
	os << "Mean" << "\t";
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << obj.mean(cat) << "\t";
	}	
	os << "\n";
	os << "Var" << "\t";
	for (unsigned int cat = 0; cat < obj.data.size(); cat++) 
	{
		os << obj.var(cat) << "\t";
	}	
	os << "\n";
	return(os);
}



////////////////////////// InvertedMStat ///////////////////////////////

// Inverted multivariate stats

InvertedMStat::InvertedMStat(const vector<vector<double> > & vec_vec_d)
	: MultivariateStat(transpose_double_matrix(vec_vec_d))
{
}

/* Exactly the same code as in Phenotype (:-@!!! )*/
vector<vector<double> > InvertedMStat::transpose_double_matrix(const vector<vector<double> > & vec_vec_d) 
{
	assert(!vec_vec_d.empty());
	
	unsigned int dim = vec_vec_d[0].size();
	vector<vector<double> > ans;
	for (unsigned int k = 0; k < dim; k++) 
	{
		vector<double> tmp;
		for (unsigned int i = 0; i < vec_vec_d.size(); i++) 
		{
			if (k==0) assert(vec_vec_d[i].size() == dim);
			tmp.push_back(vec_vec_d[i][k]);
		}
		ans.push_back(tmp);
	}
	return(ans);
}

