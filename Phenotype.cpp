// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007-2017 Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>
// Copyright 2014	   Estelle Rünneburger <estelle.runneburger@legs.cnrs-gif.fr>		

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "Phenotype.h"
#include "OutputFormat.h"

#include <cassert> // for assert

class PhenoBase;

using namespace std;

void outformat(ostream & out, const Phenotype & pheno,
                      unsigned int width, unsigned int precision, string sep)
{
    pheno.pheno_ptr->outformat(out, width, precision, sep);
}

void outformat2(ostream & out, const Phenotype & pheno,
                      unsigned int width, unsigned int precision, string sep)
{
    pheno.pheno_ptr->outformat2(out, width, precision, sep);
}

/////////////////////// Phenotype /////////////////////////////
Phenotype::Phenotype()
{
}

Phenotype::Phenotype(const Phenotype& p)
{
    if (p.pheno_ptr != nullptr)
        pheno_ptr = p.pheno_ptr->clone();
}

Phenotype::~Phenotype()
{
    // No need to delete anything because of the 
}

Phenotype& Phenotype::operator=(const Phenotype& p)
{
    if (this != &p) {
        pheno_ptr = move(p.pheno_ptr->clone());
    }
    return *this;
}

Phenotype::Phenotype(const PhenoBase & pb)
    : pheno_ptr(pb.clone())
{
}

Phenotype::Phenotype(pheno_type bp)
    : pheno_ptr(new PhenoScalar(bp))
{
}

Phenotype::Phenotype(size_t s)
    : pheno_ptr(new PhenoVector(s))
{
}

Phenotype::Phenotype(const std::vector<pheno_type>& pc)
    : pheno_ptr(new PhenoVector(pc))
{
}

Phenotype::Phenotype(const std::vector<pheno_type>& pc1, const std::vector<pheno_type>& pc2)
    : pheno_ptr(new PhenoTranscriptome(pc1, pc2))
{
}


Phenotype::Phenotype(std::unique_ptr<PhenoBase> pp)
// private constructor
    : pheno_ptr(pp->clone())
{
}

Phenotype::Phenotype(std::unique_ptr<const PhenoBase> pp)
// private constructor
    : pheno_ptr(pp->clone())
{
}

Phenotype Phenotype::sum(const vector<Phenotype>& vp)
{
    auto first = vp.begin();
    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const

    for (++first; first != vp.end(); ++first) // don't add again the first element!
        result->add(*(first->pheno_ptr)); 
        
    return result;    
}


Phenotype Phenotype::sumsq(const vector<Phenotype>& vp)
{
    auto first = vp.begin();
    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const
    result->square();
    
    for (++first; first != vp.end(); ++first) // don't add again the first element!
    {
        auto tmp = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone());
        tmp->square();
        result->add(*tmp);
     }
    return result;    
}

Phenotype Phenotype::mean(const vector<Phenotype>& vp)
{
    auto result = Phenotype::sum(vp).pheno_ptr;
    result->divide(vp.size());
    return result;
}

Phenotype Phenotype::var(const vector<Phenotype>& vp)
{
    auto mm     = Phenotype::mean(vp).pheno_ptr;
    mm->square();
    auto result = Phenotype::sumsq(vp).pheno_ptr;
    result->divide(vp.size());
    result->remove(*mm);
    return result;
}

bool Phenotype::is_defined() const
{
    return pheno_ptr != nullptr;
}

size_t Phenotype::dimensionality() const
{
    return pheno_ptr->dimensionality();
}

pheno_type Phenotype::operator[] (size_t pos) const
{
    return (*pheno_ptr)[pos];
}

pheno_type& Phenotype::operator[] (size_t pos)
{
    return (*pheno_ptr)[pos];
}

pheno_type Phenotype::get_pheno(size_t pos) const
{
    return pheno_ptr->get_pheno(pos);
}

pheno_type Phenotype::get_pheno2(size_t pos) const
{
    return pheno_ptr->get_pheno2(pos);
}
