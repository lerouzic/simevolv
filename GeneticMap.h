#ifndef GENETICMAP_H_INCLUDED
#define GENETICMAP_H_INCLUDED

#include "Parameters.h"

#include <iostream>
#include <vector>



class GeneticMap
	{
	public:
	    //constructors/destructor
	    GeneticMap();
	    GeneticMap(const ParameterSet&);
	
	    //operator overload
	    friend std::ostream& operator<< (std::ostream&, const GeneticMap&);
	
	    //functions
	    int nb_loc() const;
	    double recombination_rate(int loc1, int loc2 = -1) const;
	
	protected:
	    std::vector<double> recrate;
};


#endif // GENETICMAP_H_INCLUDED

