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





#endif
