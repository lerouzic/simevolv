#ifndef STATISTICS_H_INCLUDED
#define STATISTICS_H_INCLUDED

#include <vector>

class UnivariateStat
{
	public:
		UnivariateStat(const std::vector<double> &);
		~UnivariateStat();
		void initialize();
	
		double mean() const;
		double var() const;
	
	protected:
		const std::vector<double> data;
		
		double sum_i;
		double sum_i2;
};





#endif
