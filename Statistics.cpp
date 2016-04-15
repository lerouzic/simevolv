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
	, sum_log_i(0.0)
	, sum_log_i2(0.0)
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
			
			double log_data_i;
			if(data[i] < EXPMINLOG)
				log_data_i = MINLOG;
			else
				log_data_i = log(data[i]);
				
			sum_log_i += log_data_i;
			sum_log_i2 += log_data_i*log_data_i;
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

double UnivariateStat::mean_log() const {
	return(sum_log_i / static_cast<double>(data.size()));
}

double UnivariateStat::var_log() const {
	double mm_log = mean_log();
	return(sum_log_i2/static_cast<double>(data.size()) - mm_log*mm_log);
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
		sum_log_ij.push_back(tmp_sumij);
	}
	
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		
		// vector of sum_i
		double tmp_sum = 0.0;
		double tmp_sum_log = 0.0;
		for (unsigned int k = 0; k < data[i].size(); k++) 
		{
			tmp_sum += data[i][k];
			if (data[i][k] < EXPMINLOG)
				tmp_sum_log += MINLOG;
			else
				tmp_sum_log += log(data[i][k]);
		}
		sum_i.push_back(tmp_sum);
		sum_log_i.push_back(tmp_sum_log);
		
		// matrix of sum_ij
		for (unsigned int j = i; j < data.size(); j++) 
		{
			double tmp_sumij = 0.0;
			double tmp_log_sumij = 0.0;
			for (unsigned int k = 0; k < data[i].size(); k++)
			{
				tmp_sumij += data[i][k] * data[j][k];
				double log_data_ik = MINLOG;
				double log_data_jk = MINLOG;
				if (data[i][k] > EXPMINLOG)
					log_data_ik = log(data[i][k]);
				if (data[j][k] > EXPMINLOG)
					log_data_jk = log(data[j][k]);
				tmp_log_sumij += log_data_ik * log_data_jk;
			}
			sum_ij[i][j] = tmp_sumij;
			sum_log_ij[i][j] = tmp_log_sumij;
			if (i != j) { // does not harm, but useless when i == j
				sum_ij[j][i] = tmp_sumij; 
				sum_log_ij[j][i] = tmp_log_sumij;
			}
		}
	}
}

vector<double> MultivariateStat::means() const
{
	vector<double> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = mean(i);
	}
	return(ans);
}

vector<double> MultivariateStat::means_log() const
{
	vector<double> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = mean_log(i);
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

vector<double> MultivariateStat::vars_log() const
{
	vector<double> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = var_log(i);
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

vector<vector<double> > MultivariateStat::vcov_log() const
{
	vector<vector<double> > ans (data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		vector<double> tmp (data.size());
		for (unsigned int j = 0; j < data.size(); j++) 
		{
			tmp[j] = cov_log(i, j);
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

double MultivariateStat::mean_log(unsigned int i) const
{
	assert (i < data.size());
	return(sum_log_i[i] / static_cast<double>(data[i].size()));
}

double MultivariateStat::var(unsigned int i) const 
{
	assert (i < data.size());
	return(sum_ij[i][i]/static_cast<double>(data[i].size()) - mean(i)*mean(i));
}

double MultivariateStat::var_log(unsigned int i) const 
{
	assert (i < data.size());
	return(sum_log_ij[i][i]/static_cast<double>(data[i].size()) - mean_log(i)*mean_log(i));
}

double MultivariateStat::cov(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(sum_ij[i][j]/static_cast<double>(data[i].size()) - mean(i)*mean(j));
}

double MultivariateStat::cov_log(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(sum_log_ij[i][j]/static_cast<double>(data[i].size()) - mean_log(i)*mean_log(j));
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

////////////////////////// FastIMStat //////////////////////////////////

// Fast (and non-exhaustive) version if InvertedMStat, just for 
// internal calculation in the Network models

FastIMStat::FastIMStat(const vector<vector<double>> & dd)
	: timesteps(dd.size())
{
	// data is a vector of expression values. 
	assert(!dd.empty());
	assert(!dd[0].empty());
	size_t numphen = dd[0].size();
	for (size_t i = 0; i < numphen; i++) {
		double sx = 0.0;
		double sx2 = 0.0;
		for (size_t j = 0; j < timesteps; j++) {
			sx += dd[j][i];
			sx2 += dd[j][i]*dd[j][i];
		}
		sumx.push_back(sx);
		sumx2.push_back(sx2);
	}
}

vector<double> FastIMStat::means() const 
{
	vector<double> ans;
	for (size_t j = 0; j < sumx.size(); j++) {
		ans.push_back(sumx[j]/timesteps);
	}
	return(ans);
}

vector<double> FastIMStat::vars() const
{
	vector<double> ans;
	for (size_t j = 0; j < sumx.size(); j++) {
		double m = sumx[j]/timesteps;
		ans.push_back(sumx2[j]/timesteps - m*m);
	}
	return(ans);
}
