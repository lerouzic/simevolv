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


/* Univariate Statistics: provides the mean and variance of a vector of double */
class UnivariateStat
{
	public:
		// constructors/destructors
		UnivariateStat(const std::vector<double> &);
		
		// get results
		double mean() const;
		double var() const;
		
	protected:
		// initialization
		void initialize();
		
		// data and 
		std::vector<double> data;
		double sum_i;
		double sum_i2;
};


/* Multivariate statistics. The input is a matrix of doubles (provided as a vector of vectors)
   The function provides: the mean and variance for each entry, as well as pairwise covariances, 
   correlations, and the slope of regressions between pairs of variables */
class MultivariateStat
{
	public:
		// constructors/destructors
		MultivariateStat(const std::vector<std::vector<double> > &);

		// get infos: dimensions
		unsigned int dim1() const {return(data.size());}
		unsigned int dim2() const {if (dim1() > 0) return(data[0].size()); else return(0);}

		// get results: vectorized
		std::vector<double> means() const;
		std::vector<double> vars() const;
		std::vector<std::vector<double> > vcov() const;
		
		// ... or as scalars by providing the indexes.
		double mean(unsigned int) const;
		double var(unsigned int) const;
		double cov(unsigned int, unsigned int) const;
		double cor(unsigned int, unsigned int) const;
		double r2(unsigned int, unsigned int) const;
		double regression_slope(unsigned int, unsigned int) const; //param 1 = a*param 2 + b
		
		// output
		friend std::ostream & operator << (std::ostream &, const MultivariateStat &);
		
	protected:
		// initialization
		void initialize();
	
		const std::vector<std::vector<double> > data;
		std::vector<std::vector<double> > sum_ij;
		std::vector<double> sum_i;
};


class InvertedMStat: public MultivariateStat
{
	public:
		InvertedMStat(const std::vector<std::vector<double> > &);
	
	protected:
		// functions
		static std::vector<std::vector<double> > transpose_double_matrix(const std::vector<std::vector<double> > &);
};

class FastIMStat
{
	public: 
		FastIMStat(const std::vector<std::vector<double>> &);
		
		std::vector<double> means() const;
		std::vector<double> vars() const;
		
	protected:
		size_t timesteps;
		std::vector<double> sumx;
		std::vector<double> sumx2;
};

#endif // STATISTICS_H_INCLUDED
