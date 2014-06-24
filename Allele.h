// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef ALLELE_H_INCLUDED
#define ALLELE_H_INCLUDED

#include <vector>
#include <memory>

#include "Parameters.h"



class Allele
{
    friend class Haplotype;
    friend class Architecture; 
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiMasel;
    friend class ArchiSiegal; 
    friend class ArchiM2;
	
	public :
	    //constructors/destructor
	    Allele(const std::vector<double>);
	    Allele(const Allele &);
	    virtual ~Allele() { }
	
	    //operator overload
	    int operator== (const Allele&) const;
	    int operator!= (const Allele&) const;
	
	    //functions
	    unsigned int all_size() const;
		static std::vector<double> combine_add(const Allele &, const Allele &);
		virtual std::shared_ptr<Allele> make_mutant(double mutsd) const;
	
	protected :
	    std::vector<double> allele;
};

class Allele_zero: public Allele
{
	friend class Haplotype;
    friend class Architecture; 
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiMasel;
    friend class ArchiSiegal; 
    friend class ArchiM2;
    
	public:
	// This is a 'normal Allele' but sites with 0.0 values cannot mutate
	Allele_zero(const std::vector<double>);
	Allele_zero(const Allele_zero &);
	virtual ~Allele_zero() { }
	
	virtual std::shared_ptr<Allele> make_mutant(double mutsd) const;
};


#endif // ALLELE_H_INCLUDED
