#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

#include <vector>


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
		
	protected:
		const std::vector<std::vector<double> > data;
		std::vector<std::vector<double> > sum_ij;
		std::vector<double> sum_i;
};



#endif
