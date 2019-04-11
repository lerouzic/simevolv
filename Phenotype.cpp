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
#include "Parconst.h"
#include "OutputFormat.h"

#include <cassert> // for assert

class PhenoBase;

using namespace std;


/////////////////////// Phenotype /////////////////////////////
Phenotype::Phenotype()
{
}

Phenotype::Phenotype(const Phenotype& p)
{
    if (p.is_defined())
        pheno_ptr = p.pheno_ptr->clone();
}

Phenotype::~Phenotype()
{
    // No need to delete anything because of the unique_ptr
}

Phenotype& Phenotype::operator=(const Phenotype& p)
{
    if (p.is_defined() && (this != &p)) {
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
    assert(first->is_defined());
    
    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const

    for (++first; first != vp.end(); ++first) // don't add again the first element!
    {
        assert(first->is_defined());
        result->add(*(first->pheno_ptr)); 
    }
        
    return result;    
}


Phenotype Phenotype::sumsq(const vector<Phenotype>& vp)
{
    auto first = vp.begin();
    assert(first->is_defined());
     
    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const
    result->square();
    
    for (++first; first != vp.end(); ++first) // don't add again the first element!
    {
        assert(first->is_defined());
        auto tmp = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone());
        tmp->square();
        result->add(*tmp);
     }
    return result;    
}

Phenotype Phenotype::sumprodi(const vector<Phenotype>& vp, size_t i)
{
    auto first = vp.begin();
    assert(first->is_defined());

    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const
    result->multiply_by_index(i);

    for (++first; first != vp.end(); ++first) // don't add again the first element!
    {
        assert(first->is_defined());
        auto tmp = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone());
        tmp->multiply_by_index(i);
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

Phenotype Phenotype::vcov(const vector<Phenotype>& vp, size_t i)
{
    // i is the index of the phenotypic item (trait) against which (co)variances should be calculated. 
    auto mm = Phenotype::mean(vp).pheno_ptr;
    mm->multiply_by_index(i); // each element in mm is multiplied by mm[i]

    auto result = Phenotype::sumprodi(vp, i).pheno_ptr;
    result->divide(vp.size());
    result->remove(*mm);
    return(result);
}

bool Phenotype::is_defined() const
{
    return pheno_ptr != nullptr;
}

size_t Phenotype::dimensionality() const
{
    if(!is_defined()) return 0;
    
    return pheno_ptr->dimensionality();
}

pheno_type Phenotype::operator[] (size_t pos) const
{
    assert(is_defined());
    return (*pheno_ptr)[pos];
}

pheno_type& Phenotype::operator[] (size_t pos)
{
    assert(is_defined());
    return (*pheno_ptr)[pos];
}

pheno_type Phenotype::get_pheno(size_t pos) const
{
    assert(is_defined());
    return pheno_ptr->get_pheno(pos);
}

pheno_type Phenotype::get_pheno2(size_t pos) const
{
    assert(is_defined());
    return pheno_ptr->get_pheno2(pos);
}

void Phenotype::scale_transform(const string & option)
{
	if (option == ST_none) {
		// nothing to do
	} else if (option == ST_log) { // log transformation
		pheno_ptr->log_transform();
	} else if (option == ST_logit) { // logit(p) = log(p/(1-p))
		pheno_ptr->logit_transform();
	} else if (option == ST_m1101) { // from [-1; 1] to [0; 1]
		pheno_ptr->m1101_transform();
	} else if (option == ST_invlogit) { // invlogit(x) = 1/(1+exp(-x))
		pheno_ptr->invlogit_transform();
	} else {
		cerr << "Scale transformation error: option " << option << " unknown." << endl;
		exit(EXIT_FAILURE);
	}
}
