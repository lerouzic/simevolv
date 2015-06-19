// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>
// Copyright 2015      Christine Mayer     <christine.mayer@ibv.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef ARCHIBoolean_H_INCLUDED
#define ARCHIBoolean_H_INCLUDED

#include "Architecture.h"

#include <iostream>

class ArchiBoolean : public Architecture
{
    public :
    //constructors/destructor
    ArchiBoolean() = delete;
    ArchiBoolean(const Architecture&) = delete;
    ArchiBoolean(const ParameterSet&);
    virtual ~ArchiBoolean() {}
    
    // operator overload
    friend std::ostream& operator << (std::ostream&, const Architecture&);
    
    // getters
    virtual double get_bucket_matrix(unsigned int, unsigned int) const;
    virtual double get_logic_operator(unsigned int) const;
    
    // setters
    virtual void set_bucket_matrix(unsigned int, std::vector<double>);
    
    // debug/check
   // virtual std::string print_epsilon2() const;
    
    virtual Phenotype phenotypic_value(const Genotype&) const;
    virtual std::shared_ptr<Allele> allele_init(const ParameterSet &, unsigned int loc = 0) const;
    virtual std::shared_ptr<Allele> allele_mutation(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
    virtual std::shared_ptr<Allele> allele_mutation_test(const std::shared_ptr<Allele>, unsigned int loc = 0) const;
    
    protected :
    std::vector<std::vector<double>> bucket_matrix;
    std::vector<double> logic_operator;
    unsigned int ploc;
    std::string type;
};

#endif // ARCHIBoolean_H_INCLUDED