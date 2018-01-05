// Copyright 2014       Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

#include <vector>
#include <iostream>


const double MINLOG = -20.0;
constexpr double EXPMINLOG = 2.061154e-9; // change both at the same time!

/* Univariate Statistics: provides the mean and variance of a vector of double */
template<typename T>
class UnivariateStat
{
	public:
		// constructors/destructors
		UnivariateStat(const std::vector<T> &);
		
		// get results
		T mean() const;
		T var() const;
		
		T mean_log() const;
		T var_log() const;
		
	protected:
		// initialization
		void initialize();
		
		// data and buffers
		std::vector<T> data;
		T sum_i;
		T sum_i2;
		T sum_log_i;
		T sum_log_i2;
};


/* Multivariate statistics. The input is a matrix of doubles (provided as a vector of vectors)
   The function provides: the mean and variance for each entry, as well as pairwise covariances, 
   correlations, and the slope of regressions between pairs of variables */
template<typename T>
class MultivariateStat
{
	public:
		// constructors/destructors
		MultivariateStat(const std::vector<std::vector<T>> &);

		// get infos: dimensions
		unsigned int dim1() const {return(data.size());}
		unsigned int dim2() const {if (dim1() > 0) return(data[0].size()); else return(0);}

		// get results: vectorized
		std::vector<T> means() const;
		std::vector<T> means_log() const;		
		std::vector<T> vars() const;
		std::vector<T> vars_log() const;		
		std::vector<std::vector<T> > vcov() const;		
		std::vector<std::vector<T>> vcov_log() const;
		
		// ... or as scalars by providing the indexes.
		T mean(unsigned int) const;
		T mean_log(unsigned int) const;
		T var(unsigned int) const;
		T var_log(unsigned int) const;
		T cov(unsigned int, unsigned int) const;
		T cov_log(unsigned int, unsigned int) const;
		
		T cor(unsigned int, unsigned int) const;
		T r2(unsigned int, unsigned int) const;
		T regression_slope(unsigned int, unsigned int) const; //param 1 = a*param 2 + b
		
		// output
        template<typename T2>
		friend std::ostream & operator << (std::ostream &, const MultivariateStat<T2> &);
		
	protected:
		// initialization
		void initialize();
	
		const std::vector<std::vector<T>> data;
		std::vector<std::vector<T>> sum_ij;
		std::vector<std::vector<T>> sum_log_ij;
		std::vector<T> sum_i;
		std::vector<T> sum_log_i;
};

template<typename T>
class InvertedMStat: public MultivariateStat<T>
{
	public:
		InvertedMStat(const std::vector<std::vector<T>> &);
	
	protected:
		// functions
		static std::vector<std::vector<T>> transpose_double_matrix(const std::vector<std::vector<T>> &);
};

template<typename T>
class FastIMStat
{
	public: 
		FastIMStat(const std::vector<std::vector<T>> &);
		
		std::vector<T> means() const;
		std::vector<T> vars() const;
		
	protected:
		size_t timesteps;
		std::vector<T> sumx;
		std::vector<T> sumx2;
};

#include "Statistics.cpp" // Template class initialization

#endif // STATISTICS_H_INCLUDED
