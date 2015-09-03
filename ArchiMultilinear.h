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



#ifndef ARCHIMULTILINEAR_H_INCLUDED
#define ARCHIMULTILINEAR_H_INCLUDED

#include "Architecture.h"

#include <iostream>

class ArchiMultilinear : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiMultilinear() = delete;
	    ArchiMultilinear(const Architecture&) = delete;
	    ArchiMultilinear(const ParameterSet&);
	    virtual ~ArchiMultilinear() {}
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    // getters
	    virtual double get_epsilon2(unsigned int, unsigned int) const;
	    virtual double get_epsilon3(unsigned int, unsigned int, unsigned int) const;
	    
	    // setters
	    virtual void set_epsilon2(unsigned int, unsigned int, double);
	    virtual void set_epsilon3(unsigned int, unsigned int, unsigned int, double);
	    
	    // debug/check
	    virtual std::string print_epsilon2() const;
	    virtual std::string print_epsilon3() const;
	
	    virtual Phenotype phenotypic_value(const Genotype&, bool envir) const;
	
	protected :
	    std::vector<std::vector<double>> epsilon2;
	    std::vector<std::vector<std::vector<double>>> epsilon3;
	    bool flag_epistasis2;
	    bool flag_epistasis3;
	    
		// flags to speed up calculations
	    virtual bool is_epistasis() const {return((is_epistasis2()) || (is_epistasis3()));}
	    virtual bool is_epistasis2() const {return(flag_epistasis2);}
	    virtual bool is_epistasis3() const {return(flag_epistasis3);}	    
};

#endif // ARCHIMULTILINEAR_H_INCLUDED
