#include <cassert>
#include <iostream>

#include "Statistics.h"

using namespace std;

UnivariateStat::UnivariateStat(const vector<double> & vv)
	: data(vv)
{
	if (data.size() < 2) {
		cerr << "A data set of size " << data.size() << " is not eligible for statistics." << endl;
		assert("Terminate.");
	}
	initialize();
}

UnivariateStat::~UnivariateStat()
{
	
}

void UnivariateStat::initialize()
{
	for (unsigned int i = 0; i < data.size(); i++) {
		sum_i += data[i];
		sum_i2 += data[i]*data[i];
	}
}

double UnivariateStat::mean() const {
	return(sum_i / static_cast<double>(data.size()));
}

double UnivariateStat::var() const {
	double mm = mean();
	return(sum_i2/static_cast<double>(data.size()) - mm*mm);
}

