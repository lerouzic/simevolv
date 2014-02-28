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



#ifndef ARCHIMULTILINEAR_H_INCLUDED
#define ARCHIMULTILINEAR_H_INCLUDED

#include "Architecture.h"



class ArchiMultilinear : public Architecture
{
	public :
	    //constructors/destructor
	    ArchiMultilinear();
	    ArchiMultilinear(const Architecture&);
	    ArchiMultilinear(const ParameterSet&);
	    ~ArchiMultilinear() {}
	
	    // operator overload
	    friend std::ostream& operator << (std::ostream&, const Architecture&);
	
	    //functions
	    double get_epsilon2(unsigned int, unsigned int) const;
	    double get_epsilon3(unsigned int, unsigned int, unsigned int) const;
	    void set_epsilon2(unsigned int, unsigned int, double);
	    void set_epsilon3(unsigned int, unsigned int, unsigned int, double);
	    std::string print_epsilon2() const;
	    std::string print_epsilon3() const;
	
	    bool is_epistasis() const {return((is_epistasis2()) || (is_epistasis3()));}
	    bool is_epistasis2() const {return(flag_epistasis2);}
	    bool is_epistasis3() const {return(flag_epistasis3);}
	
	    Phenotype phenotypic_value(const Genotype&) const;
	
	protected :
	    std::vector<std::vector<double> > epsilon2;
	    std::vector<std::vector<std::vector<double> > > epsilon3;
	    bool flag_epistasis2;
	    bool flag_epistasis3;

};


#endif // ARCHIMULTILINEAR_H_INCLUDED
