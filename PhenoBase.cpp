// Copyright 2017 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/




#include "PhenoBase.h"
#include "OutputFormat.h"

#include <cassert> // for assert

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/export.hpp>
BOOST_SERIALIZATION_ASSUME_ABSTRACT(PhenoBase)
BOOST_CLASS_EXPORT(PhenoScalar)
BOOST_CLASS_EXPORT(PhenoVector)
BOOST_CLASS_EXPORT(PhenoTranscriptome)
#endif

using namespace std;

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

void PhenoScalar::multiply_by_index(size_t i)
{
    assert(i < dimensionality());
    pheno *= pheno; // this is the same as squaring the phenotype
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
    : pheno(std::vector<pheno_type>())
{
}

PhenoVector::PhenoVector(size_t s)
    : pheno(std::vector<pheno_type>(s))
{
}

PhenoVector::PhenoVector(const PhenoVector& p)
    : pheno(p.pheno)
{
}

PhenoVector::PhenoVector(const std::vector<pheno_type>& v)
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

void PhenoVector::multiply_by_index(size_t i)
{
    assert(i < dimensionality());
    //~ for (auto &p : pheno)
        //~ p *= pheno[i];
    for (unsigned int j = 0; j < pheno.size(); j++) {
        pheno[j] *= pheno[i];
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
    : expression(std::vector<pheno_type>())
    , stability(std::vector<pheno_type>())
{
}

PhenoTranscriptome::PhenoTranscriptome(size_t n)
    : expression(std::vector<pheno_type>(n))
    , stability(std::vector<pheno_type>(n))
{
}

PhenoTranscriptome::PhenoTranscriptome(const PhenoTranscriptome& p)
    : expression(p.expression)
    , stability(p.stability)
{
}

PhenoTranscriptome::PhenoTranscriptome(const std::vector<pheno_type>& v1, const std::vector<pheno_type>& v2)
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

void PhenoTranscriptome::multiply_by_index(size_t i)
{
    for (unsigned int j = 0; j < expression.size(); j++) {
        expression[j] *= expression[i];
    }
    for (unsigned int j = 0; j < stability.size(); j++) {
        stability[j] *= stability[i];
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
