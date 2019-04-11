// Copyright 2017 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef PHENOBASE_H_INCLUDED
#define PHENOBASE_H_INCLUDED

#include "types.h"

#include <memory>
#include <vector>
#include <iostream>
#include <string>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

class PhenoBase    // This is the base class, it is just and interface (virtual pure)
{
	public:
		virtual ~PhenoBase() { }
        
        virtual std::unique_ptr<PhenoBase> clone() const = 0;
        virtual void add(const PhenoBase &)    = 0;
        virtual void remove(const PhenoBase &) = 0;
        virtual void divide(unsigned int)      = 0;
        virtual void square()                  = 0;
        virtual void multiply_by_index(std::size_t) = 0;
        
        virtual std::size_t dimensionality() const = 0;
        virtual pheno_type operator[] (std::size_t) const  = 0;
        virtual pheno_type& operator[](std::size_t)        = 0;
        virtual pheno_type get_pheno(std::size_t) const;
        virtual pheno_type get_pheno2(std::size_t) const;      
        
        virtual void log_transform()      = 0;
        virtual void logit_transform()    = 0;
        virtual void m1101_transform()    = 0;
        virtual void invlogit_transform() = 0;
        
        virtual void outformat(std::ostream&, unsigned int, unsigned int, std::string) const = 0;
        virtual void outformat2(std::ostream&, unsigned int, unsigned int, std::string) const;
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
        }
        #endif
};

class PhenoScalar: public PhenoBase {
  // The phenotype is a single scalar
  
  // Phenotypes should be canonical classes
  public: 
    PhenoScalar();
    PhenoScalar(const PhenoScalar&);
    PhenoScalar(pheno_type);
    PhenoScalar& operator=(const PhenoScalar&);
    ~PhenoScalar() { }
    
    std::unique_ptr<PhenoBase> clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    void multiply_by_index(std::size_t);
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t);
    
    void log_transform();
    void logit_transform();
    void m1101_transform();
    void invlogit_transform();
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;    
          
  protected:
    
    PhenoScalar  operator +  (const PhenoScalar&) const;
    PhenoScalar& operator += (const PhenoScalar&);
    PhenoScalar  operator -  (const PhenoScalar&) const;
    PhenoScalar& operator -= (const PhenoScalar&);    
    PhenoScalar  operator /  (unsigned int) const;
    PhenoScalar& operator /= (unsigned int);  
  
    pheno_type pheno;
    
  private:
  #ifdef SERIALIZATION_TEXT
	friend class boost::serialization::access;
	template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<PhenoBase>(*this);  
      ar & pheno;
    }
  #endif    
};

class PhenoVector: public PhenoBase {
  // The phenotype is a vector
  
  public: 
    PhenoVector(); 
    PhenoVector(size_t);
    PhenoVector(const PhenoVector&);
    PhenoVector(const std::vector<pheno_type>&);
    PhenoVector& operator = (const PhenoVector& p);
    ~PhenoVector() { }    
    
    std::unique_ptr<PhenoBase> clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    void multiply_by_index(std::size_t);
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t);
    
    void log_transform();
    void logit_transform();
    void m1101_transform();
    void invlogit_transform();
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;      
        
  protected:
    PhenoVector  operator +  (const PhenoVector&) const;
    PhenoVector& operator += (const PhenoVector&);
    PhenoVector  operator -  (const PhenoVector&) const;
    PhenoVector& operator -= (const PhenoVector&);    
    PhenoVector  operator /  (unsigned int) const;
    PhenoVector& operator /= (unsigned int);  
  
    std::vector<pheno_type> pheno;
    
  private:
  #ifdef SERIALIZATION_TEXT
	friend class boost::serialization::access; 
	template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<PhenoBase>(*this);         
      ar & pheno;
    }
  #endif     
};

class PhenoTranscriptome: public PhenoBase {
  // The phenotype is a transcriptome
  //   - gene expression levels
  //   - gene expression stability
  
  public: 
    PhenoTranscriptome();
    PhenoTranscriptome(size_t);
    PhenoTranscriptome(const PhenoTranscriptome&);
    PhenoTranscriptome(const std::vector<pheno_type>&, const std::vector<pheno_type>&);
    PhenoTranscriptome& operator= (const PhenoTranscriptome&);
    ~PhenoTranscriptome() { }
    
    std::unique_ptr<PhenoBase> clone() const;     
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    void multiply_by_index(std::size_t);
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t); 
    pheno_type  get_pheno2(std::size_t) const;
    
    void log_transform();
    void logit_transform();
    void m1101_transform();
    void invlogit_transform();
    
    void outformat (std::ostream&, unsigned int, unsigned int, std::string) const; 
    void outformat2(std::ostream&, unsigned int, unsigned int, std::string) const;             
    
  protected: 
     
    PhenoTranscriptome  operator +  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator += (const PhenoTranscriptome&); 
    PhenoTranscriptome  operator -  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator -= (const PhenoTranscriptome&);     
    PhenoTranscriptome  operator /  (unsigned int) const;
    PhenoTranscriptome& operator /= (unsigned int); 
  
    std::vector<pheno_type> expression;
    std::vector<pheno_type> stability;
    
  private:
  #ifdef SERIALIZATION_TEXT
	friend class boost::serialization::access;
	template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<PhenoBase>(*this);  
      ar & expression;
      ar & stability;
    }
  #endif     
};

#endif // # PHENOTYPE_H_INCLUDED
