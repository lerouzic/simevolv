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

Phenotype::Phenotype(const Phenotype::basic_container& pc)
    : pheno_ptr(new PhenoVector(pc))
{
}

Phenotype::Phenotype(const Phenotype::basic_container& pc1, const Phenotype::basic_container& pc2)
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

Phenotype Phenotype::sum(const pheno_container& vp)
{
    auto first = vp.begin();
    auto result = std::unique_ptr<PhenoBase>(first->pheno_ptr->clone()); // because *first is const

    for (++first; first != vp.end(); ++first) // don't add again the first element!
        result->add(*(first->pheno_ptr)); 
        
    return result;    
}


Phenotype Phenotype::sumsq(const pheno_container& vp)
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

Phenotype Phenotype::mean(const pheno_container& vp)
{
    auto result = Phenotype::sum(vp).pheno_ptr;
    result->divide(vp.size());
    return result;
}

Phenotype Phenotype::var(const pheno_container& vp)
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

/////////////////////// PhenoBase   ///////////////////////////

pheno_type PhenoBase::get_pheno(size_t pos) const
{
    return (*this)[pos];
}

pheno_type PhenoBase::get_pheno2(size_t pos) const
{
    assert("This should never happen: calling get_pheno2 on an irrelevant phenotype.\n");
    return 0.0;
}

void PhenoBase::outformat2(ostream & out, unsigned int width, unsigned int precision, string sep) const
{
    assert("This should never happen: calling outformat2 on an irrelevant phenotype.\n");
}    
/////////////////////// PhenoScalar ///////////////////////////

PhenoScalar::PhenoScalar() 
    : pheno(0.0)
{
}

PhenoScalar::PhenoScalar(const PhenoScalar& p)
    : pheno(p.pheno)
{
}

PhenoScalar::PhenoScalar(pheno_type v)
    : pheno(v)
{
}

PhenoScalar& PhenoScalar::operator = (const PhenoScalar& p)
{ 
    if (this != &p) {
        pheno = p.pheno;
    }
    return *this;
}

std::unique_ptr<PhenoBase> PhenoScalar::clone() const
{
    return std::unique_ptr<PhenoBase>(new PhenoScalar(*this));
}

PhenoScalar& PhenoScalar::operator += (const PhenoScalar& p)
{
    pheno += p.pheno;
    return *this;
}

PhenoScalar PhenoScalar::operator + (const PhenoScalar& p) const
{
    return PhenoScalar(*this) += p;
}

PhenoScalar& PhenoScalar::operator -= (const PhenoScalar& p)
{
    pheno -= p.pheno;
    return *this;
}

PhenoScalar PhenoScalar::operator - (const PhenoScalar& p) const
{
    return PhenoScalar(*this) -= p;
}


PhenoScalar& PhenoScalar::operator /= (unsigned int n)
{
    pheno /= n;
    return *this;
}

PhenoScalar PhenoScalar::operator / (unsigned int n) const
{
    return PhenoScalar(*this) /= n;
}

void PhenoScalar::add(const PhenoBase& p)
{
    const PhenoScalar& pp = dynamic_cast<const PhenoScalar&>(p);
    *this += pp;
}

void PhenoScalar::remove(const PhenoBase& p)
{
    const PhenoScalar& pp = dynamic_cast<const PhenoScalar&>(p);
    *this -= pp;
}

void PhenoScalar::divide(unsigned int n)
{
    *this /= n;
}

void PhenoScalar::square()
{
    pheno *= pheno;
}

size_t PhenoScalar::dimensionality() const
{
    return 1;
}

pheno_type PhenoScalar::operator[] (size_t pos) const
{
    return pheno;
}

pheno_type& PhenoScalar::operator[] (size_t pos)
{
    return pheno;
}

void PhenoScalar::outformat(ostream& out, unsigned int width, unsigned int precision, string sep) const
{
    ::outformat(out, pheno, width, precision, sep);
}

////////////////////// PhenoVector /////////////////////////////

PhenoVector::PhenoVector()
    : pheno(Phenotype::basic_container())
{
}

PhenoVector::PhenoVector(size_t s)
    : pheno(Phenotype::basic_container(s))
{
}

PhenoVector::PhenoVector(const PhenoVector& p)
    : pheno(p.pheno)
{
}

PhenoVector::PhenoVector(const Phenotype::basic_container& v)
    : pheno(v)
{
}

PhenoVector& PhenoVector::operator = (const PhenoVector & p)
{
    if (this != &p) {
        pheno = p.pheno;
    }
    return *this;
}

std::unique_ptr<PhenoBase> PhenoVector::clone() const
{
    return std::unique_ptr<PhenoBase>(new PhenoVector(*this));
}

PhenoVector& PhenoVector::operator += (const PhenoVector & p)
{
    assert(pheno.size() == p.pheno.size());
    for (unsigned int i = 0; i < pheno.size(); i++) {
        pheno[i] += p.pheno[i];
    }
    return *this;
}

PhenoVector PhenoVector::operator + (const PhenoVector & p) const
{
    return PhenoVector(*this) += p;
}

PhenoVector& PhenoVector::operator -= (const PhenoVector & p)
{
    assert(pheno.size() == p.pheno.size());
    for (unsigned int i = 0; i < pheno.size(); i++) {
        pheno[i] -= p.pheno[i];
    }
    return *this;
}

PhenoVector PhenoVector::operator - (const PhenoVector & p) const
{
    return PhenoVector(*this) -= p;
}

PhenoVector& PhenoVector::operator /= (unsigned int n)
{
    for (auto &i: pheno) {
        i /= n;
    }
    return *this;
}

PhenoVector PhenoVector::operator / (unsigned int n) const
{
    return PhenoVector(*this) /= n;
}

void PhenoVector::add(const PhenoBase& p)
{
    const PhenoVector& pp = dynamic_cast<const PhenoVector&>(p);
    *this += pp;
}

void PhenoVector::remove(const PhenoBase& p)
{
    const PhenoVector& pp = dynamic_cast<const PhenoVector&>(p);
   *this -= pp;
}

void PhenoVector::divide(unsigned int n)
{
    *this /= n;
}

void PhenoVector::square()
{
    for (unsigned int i = 0; i < pheno.size(); i++) {
        pheno[i] *= pheno[i];
    }    
}

size_t PhenoVector::dimensionality() const
{
    return pheno.size();
}
pheno_type PhenoVector::operator[] (size_t pos) const
{
    return pheno[pos];
}

pheno_type& PhenoVector::operator[] (size_t pos)
{
    return pheno[pos];
}

void PhenoVector::outformat(ostream& out, unsigned int width, unsigned int precision, string sep) const
{
    for (auto p : pheno)
        ::outformat(out, p, width, precision, sep);
}


////////////////////// class PhenoTranscriptome ///////////////////////

PhenoTranscriptome::PhenoTranscriptome() 
    : expression(Phenotype::basic_container())
    , stability(Phenotype::basic_container())
{
}

PhenoTranscriptome::PhenoTranscriptome(size_t n)
    : expression(Phenotype::basic_container(n))
    , stability(Phenotype::basic_container(n))
{
}

PhenoTranscriptome::PhenoTranscriptome(const PhenoTranscriptome& p)
    : expression(p.expression)
    , stability(p.stability)
{
}

PhenoTranscriptome::PhenoTranscriptome(const Phenotype::basic_container& v1, const Phenotype::basic_container& v2)
    : expression(v1)
    , stability(v2)
{
}

PhenoTranscriptome& PhenoTranscriptome::operator = (const PhenoTranscriptome& p)
{
    if (this != &p) {
        expression = p.expression;
        stability  = p.stability;
    }
    return *this;
}

std::unique_ptr<PhenoBase> PhenoTranscriptome::clone() const
{
    return std::unique_ptr<PhenoBase>(new PhenoTranscriptome(*this));
}

PhenoTranscriptome& PhenoTranscriptome::operator += (const PhenoTranscriptome &p)
{
    assert(expression.size() == p.expression.size());
    assert(stability.size() == p.stability.size());
    for (unsigned int i = 0; i < expression.size(); i++)
        expression[i] += p.expression[i];
    for (unsigned int i = 0; i < stability.size(); i++)
        stability[i] += p.stability[i];
    return *this;
}

PhenoTranscriptome PhenoTranscriptome::operator + (const PhenoTranscriptome &p) const
{
    return PhenoTranscriptome(*this) += p;
}

PhenoTranscriptome& PhenoTranscriptome::operator -= (const PhenoTranscriptome &p)
{
    assert(expression.size() == p.expression.size());
    assert(stability.size() == p.stability.size());
    for (unsigned int i = 0; i < expression.size(); i++)
        expression[i] -= p.expression[i];
    for (unsigned int i = 0; i < stability.size(); i++)
        stability[i] -= p.stability[i];
    return *this;
}

PhenoTranscriptome PhenoTranscriptome::operator - (const PhenoTranscriptome &p) const
{
    return PhenoTranscriptome(*this) -= p;
}

PhenoTranscriptome & PhenoTranscriptome::operator /= (unsigned int n)
{
    for (unsigned int i = 0; i < expression.size(); i++)
        expression[i] /= n;
    for (unsigned int i = 0; i < stability.size(); i++)
        stability[i] /= n;
    return *this;
}

PhenoTranscriptome PhenoTranscriptome::operator / (unsigned int n) const
{
    return PhenoTranscriptome(*this) /= n;
}

void PhenoTranscriptome::add(const PhenoBase& p)
{
    const PhenoTranscriptome& pp = dynamic_cast<const PhenoTranscriptome&>(p);
    *this += pp;
}

void PhenoTranscriptome::remove(const PhenoBase& p)
{
    const PhenoTranscriptome& pp = dynamic_cast<const PhenoTranscriptome&>(p);
    *this -= pp;
}

void PhenoTranscriptome::divide(unsigned int n)
{
    *this /= n;
}

void PhenoTranscriptome::square()
{
    for (unsigned int i = 0; i < expression.size(); i++) {
        expression[i] *= expression[i];
    }    
    for (unsigned int i = 0; i < stability.size(); i++) {
        stability[i] *= stability[i];
    }      
}

size_t PhenoTranscriptome::dimensionality() const
{
    return expression.size();
}

pheno_type PhenoTranscriptome::operator[] (size_t pos) const
{
    return expression[pos];
}

pheno_type& PhenoTranscriptome::operator[] (size_t pos)
{
    return expression[pos];
}

pheno_type PhenoTranscriptome::get_pheno2(size_t pos) const
{
    return stability[pos];
}

void PhenoTranscriptome::outformat(ostream& out, unsigned int width, unsigned int precision, string sep) const
{
    for (auto p : expression)
        ::outformat(out, p, width, precision, sep);
}

void PhenoTranscriptome::outformat2(ostream& out, unsigned int width, unsigned int precision, string sep) const
{
    for (auto p : stability)
        ::outformat(out, p, width, precision, sep);
}
