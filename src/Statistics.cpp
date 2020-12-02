// Copyright 2014       Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



// #include "Statistics.h"

#include <cassert>
#include <iostream>
#include <cmath>

// This is a .hpp file, to be compiled along with Statistics.h


////////////////////////////////// UNIVARIATE ////////////////////////////////////////

// Constructors and initialization 

template<typename T>
UnivariateStat<T>::UnivariateStat(const std::vector<T> & vv)
	: data(vv)
    , N_finite(0)
    , N_positive(0)
	, sum_i(0.0)
	, sum_i2(0.0)
	, sum_log_i(0.0)
	, sum_log_i2(0.0)
{
	if (data.size() < 2) 
	{
		// std::cerr << "A data set of size " << data.size() << " may not lead to meaningful statistical estimates." << std::endl;
        // Actually, everything behaves more or less properly when n=1
        // Variances are sample variance (sum sq /n) and can thus be computed
	}
	initialize();
}

/* computes and stores intermediate results. */
template<typename T>
void UnivariateStat<T>::initialize()
{
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		if (std::isfinite(data[i])) 
		{
            N_finite++;
			sum_i += data[i];
			sum_i2 += data[i]*data[i];
			
            if (data[i] > 0) {
                N_positive++;
                T log_data_i;
                if(data[i] < EXPMINLOG)
                    log_data_i = MINLOG;
                else
                    log_data_i = log(data[i]);
				
                sum_log_i += log_data_i;
                sum_log_i2 += log_data_i*log_data_i;
            }
		} 
	}
}


// functions

template<typename T>
T UnivariateStat<T>::mean() const {
	return(sum_i / static_cast<T>(N_finite));
}

template<typename T>
T UnivariateStat<T>::var() const {
	double mm = mean();
	return(sum_i2/static_cast<T>(N_finite) - mm*mm);
}

template<typename T>
T UnivariateStat<T>::mean_log() const {
	return(sum_log_i / static_cast<T>(N_positive));
}

template<typename T>
T UnivariateStat<T>::var_log() const {
	T mm_log = mean_log();
	return(sum_log_i2/static_cast<T>(N_positive) - mm_log*mm_log);
}


////////////////////////////////// MULTIVARIATE /////////////////////////////////

template<typename T>
MultivariateStat<T>::MultivariateStat(const std::vector<std::vector<T> > & v)
	: data(v)
    , sum_ij(std::vector<std::vector<T>>(v.size(), std::vector<T>(v.size(), 0.0)))
    , sum_log_ij(std::vector<std::vector<T>>(v.size(), std::vector<T>(v.size(), 0.0)))
    , sum_i(std::vector<T>(v.size()))
    , sum_log_i(std::vector<T>(v.size()))
{ 
	// It is a bit complicated to assert here whether the vector<vector<double> > is OK. 
	// Temporary storage variables as well as checks on the dimensions will occur afterwards, 
	// in initialize(), which computes the sums of x_i and x_i^2. 
	initialize();
}

template<typename T>
void MultivariateStat<T>::initialize()
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
	
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		
		// vector of sum_i
		for (unsigned int k = 0; k < data[i].size(); k++) 
		{
            T d_ik =  data[i][k];
			sum_i[i] += d_ik;
			if (d_ik < EXPMINLOG)
				sum_log_i[i] += MINLOG;
			else
				sum_log_i[i] += log(d_ik);
		}
		
		// matrix of sum_ij
		for (unsigned int j = i; j < data.size(); j++) 
		{
            T tmp_sum_ij = 0.0;
            T tmp_sum_log_ij = 0.0;
			for (unsigned int k = 0; k < data[i].size(); k++)
			{
                T d_ik = data[i][k];
                T d_jk = data[j][k];
				tmp_sum_ij += d_ik * d_jk;
				T log_data_ik = MINLOG;
				T log_data_jk = MINLOG;
				if (d_ik > EXPMINLOG)
					log_data_ik = log(d_ik);
				if (d_jk > EXPMINLOG)
					log_data_jk = log(d_jk);
				tmp_sum_log_ij += log_data_ik * log_data_jk;
			}
            sum_ij[i][j] = tmp_sum_ij;
            sum_log_ij[i][j] = tmp_sum_log_ij;
			if (i != j) { // does not harm, but useless when i == j
				sum_ij[j][i] = tmp_sum_ij;
				sum_log_ij[j][i] = tmp_sum_log_ij;
			}
		}
	}
}

template<typename T>
std::vector<T> MultivariateStat<T>::means() const
{
	std::vector<T> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = mean(i);
	}
	return(ans);
}

template<typename T>
std::vector<T> MultivariateStat<T>::means_log() const
{
	std::vector<T> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = mean_log(i);
	}
	return(ans);
}

template<typename T>
std::vector<T> MultivariateStat<T>::vars() const
{
	std::vector<T> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = var(i);
	}
	return(ans);
}

template<typename T>
std::vector<T> MultivariateStat<T>::vars_log() const
{
	std::vector<T> ans(data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		ans[i] = var_log(i);
	}
	return(ans);
}

template<typename T>
std::vector<std::vector<T>> MultivariateStat<T>::vcov() const
{
	std::vector<std::vector<T>> ans (data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		std::vector<T> tmp (data.size());
		for (unsigned int j = 0; j < data.size(); j++) 
		{
			tmp[j] = cov(i, j);
		}
		ans.push_back(tmp);
	}
	return(ans);
}

template<typename T>
std::vector<std::vector<T>> MultivariateStat<T>::vcov_log() const
{
	std::vector<std::vector<T>> ans (data.size());
	for (unsigned int i = 0; i < data.size(); i++) 
	{
		std::vector<T> tmp (data.size());
		for (unsigned int j = 0; j < data.size(); j++) 
		{
			tmp[j] = cov_log(i, j);
		}
		ans.push_back(tmp);
	}
	return(ans);
}

template<typename T>
T MultivariateStat<T>::mean(unsigned int i) const
{
	assert (i < data.size());
	return(sum_i[i] / static_cast<T>(data[i].size()));
}

template<typename T>
T MultivariateStat<T>::mean_log(unsigned int i) const
{
	assert (i < data.size());
	return(sum_log_i[i] / static_cast<T>(data[i].size()));
}

template<typename T>
T MultivariateStat<T>::var(unsigned int i) const 
{
	assert (i < data.size());
	return(sum_ij[i][i]/static_cast<T>(data[i].size()) - mean(i)*mean(i));
}

template<typename T>
T MultivariateStat<T>::var_log(unsigned int i) const 
{
	assert (i < data.size());
	return(sum_log_ij[i][i]/static_cast<T>(data[i].size()) - mean_log(i)*mean_log(i));
}

template<typename T>
T MultivariateStat<T>::cov(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(sum_ij[i][j]/static_cast<T>(data[i].size()) - mean(i)*mean(j));
}

template<typename T>
T MultivariateStat<T>::cov_log(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(sum_log_ij[i][j]/static_cast<T>(data[i].size()) - mean_log(i)*mean_log(j));
}

template<typename T>
T MultivariateStat<T>::cor(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cov(i, j)/(sqrt(var(i)*var(j))));
}

template<typename T>
T MultivariateStat<T>::r2(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cor(i, j)*cor(i, j));	
}

template<typename T>
T MultivariateStat<T>::regression_slope(unsigned int i, unsigned int j) const
{
	assert (i < data.size());
	assert (j < data.size());
	return(cov(i, j)/var(j));
}

template<typename T>
std::ostream & operator << (std::ostream & os, const MultivariateStat<T> & obj) 
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
template<typename T>
InvertedMStat<T>::InvertedMStat(const std::vector<std::vector<T>> & vec_vec_d)
	: MultivariateStat<T>(transpose_double_matrix(vec_vec_d))
{
}

/* Exactly the same code as in Phenotype (:-@!!! )*/
template<typename T>
std::vector<std::vector<T> > InvertedMStat<T>::transpose_double_matrix(const std::vector<std::vector<T>> & vec_vec_d) 
{
	assert(!vec_vec_d.empty());
	
	unsigned int dim = vec_vec_d[0].size();
	std::vector<std::vector<T>> ans;
	for (unsigned int k = 0; k < dim; k++) 
	{
		std::vector<T> tmp;
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

template<typename T>
FastIMStat<T>::FastIMStat(const std::vector<std::vector<T>> & dd)
	: timesteps(dd.size())
{
	// data is a vector of expression values. 
	assert(!dd.empty());
	assert(!dd[0].empty());
	size_t numphen = dd[0].size();
	for (size_t i = 0; i < numphen; i++) {
		T sx = 0.0;
		T sx2 = 0.0;
		for (size_t j = 0; j < timesteps; j++) {
			sx += dd[j][i];
			sx2 += dd[j][i]*dd[j][i];
		}
		sumx.push_back(sx);
		sumx2.push_back(sx2);
	}
}

template<typename T>
std::vector<T> FastIMStat<T>::means() const 
{
	std::vector<T> ans;
	for (size_t j = 0; j < sumx.size(); j++) {
		ans.push_back(sumx[j]/timesteps);
	}
	return(ans);
}

template<typename T>
std::vector<T> FastIMStat<T>::vars() const
{
	std::vector<T> ans;
	for (size_t j = 0; j < sumx.size(); j++) {
		T m = sumx[j]/timesteps;
		ans.push_back(sumx2[j]/timesteps - m*m);
	}
	return(ans);
}
