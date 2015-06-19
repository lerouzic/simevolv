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



#ifndef ALLELE_H_INCLUDED
#define ALLELE_H_INCLUDED

#include "Parameters.h"

#include <vector>
#include <memory>


class Allele
{
    friend class Haplotype;
    friend class Architecture; 
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
    friend class ArchiRegulatoryMatrix;
    friend class ArchiWagner;
    friend class ArchiSiegal; 
    friend class ArchiM2;
    friend class ArchiBoolean;
	
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
		static std::vector<double> combine_mean(const Allele &, const Allele &);
        static std::vector<double> combine_OR(const Allele &, const Allele &);
		
		virtual std::shared_ptr<Allele> make_mutant(double mutsd) const;
        virtual std::shared_ptr<Allele> make_boolean_mutant() const;
	
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
