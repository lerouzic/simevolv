// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef GENETICMAP_H_INCLUDED
#define GENETICMAP_H_INCLUDED

#include "types.h"
#include "Parameters.h"

#include <vector>

 #ifdef SERIALIZATION_TEXT
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif

class GeneticMap 
{
	public:
	    //constructors/destructor
	    GeneticMap();
	    GeneticMap(const GeneticMap&) = delete;
	    GeneticMap(const ParameterSet&);
	    ~GeneticMap();
	
	    //functions
	    int nb_loc() const;
	    rate_type recombination_rate(unsigned int loc) const;
	
	protected:
	    std::vector<rate_type> recrate;
	    
	private:
         #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;	
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & recrate;
        }
        #endif
};


#endif // GENETICMAP_H_INCLUDED
