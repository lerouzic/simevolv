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



#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

#include <vector>
#include <iostream>



class UnivariateStat
{
	public:
		// constructors/destructors
		UnivariateStat(const std::vector<double> &);
		~UnivariateStat();
		
		//initialization
		void initialize();
		
		//functions
		double mean() const;
		double var() const;
		
	protected:
		const std::vector<double> data;
		double sum_i;
		double sum_i2;
};



class MultivariateStat
{
	public:
		// constructors/destructors
		MultivariateStat(const std::vector<std::vector<double> > &);
		~MultivariateStat();
		
		// initilization
		void initialize();
		
		// functions
		std::vector<double> means() const;
		std::vector<double> vars() const;
		std::vector<std::vector<double> > vcov() const;
		
		double mean(unsigned int) const;
		double var(unsigned int) const;
		double cov(unsigned int, unsigned int) const;
		double cor(unsigned int, unsigned int) const;
		double r2(unsigned int, unsigned int) const;
		double regression_slope(unsigned int, unsigned int) const; //param 1 = a*param 2 + b
		
		// output
		friend std::ostream & operator << (std::ostream &, const MultivariateStat &);
		
	protected:
		const std::vector<std::vector<double> > data;
		std::vector<std::vector<double> > sum_ij;
		std::vector<double> sum_i;
};

#endif
