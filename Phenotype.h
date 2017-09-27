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



#ifndef PHENOTYPE_H_INCLUDED
#define PHENOTYPE_H_INCLUDED

#include <iostream>
#include <string>  
#include <vector>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#endif

class Phenotype;

class PhenoBase;
typedef double basic_pheno;
typedef std::vector<Phenotype> pheno_container;
typedef std::vector<basic_pheno> basic_container;

// These ones should not really be used outside of Phenotype
typedef std::unique_ptr<PhenoBase> pheno_type;
typedef std::unique_ptr<const PhenoBase> const_pheno_type;

void outformat(std::ostream &, const Phenotype &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");
void outformat2(std::ostream &, const Phenotype &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");

class Phenotype
{ // Class Phenotype encapsulates a polymorphic pointer to a PhenoBase-derived class 
  // This is supposed to be the only interface with most of the software
  
	friend void outformat(std::ostream &, const Phenotype &, 
			unsigned int width, unsigned int precision, std::string sep);  
	friend void outformat2(std::ostream &, const Phenotype &, 
			unsigned int width, unsigned int precision, std::string sep);  
              
    public:
        Phenotype();
        Phenotype(const Phenotype &);
        Phenotype(const PhenoBase &);
        ~Phenotype();
        Phenotype& operator=(const Phenotype &);
        
        // Specialized constructors, call directly a subclass (good idea or not?)
        Phenotype(basic_pheno);  // PhenoScalar
        Phenotype(size_t);       // PhenoVector
        Phenotype(const basic_container&); // PhenoVector
        Phenotype(const basic_container&, const basic_container&); // PhenoTranscriptome
        
        static Phenotype sum (const pheno_container&);
        static Phenotype sumsq(const pheno_container&);
        static Phenotype mean(const pheno_container&);        
        static Phenotype var(const pheno_container&);

    // Necessary for compatibility of the old code, should probably be avoided   
        bool is_defined() const;
        std::size_t dimensionality() const;
        basic_pheno operator[] (std::size_t) const;
        basic_pheno &operator[] (std::size_t);
        virtual basic_pheno get_pheno(std::size_t) const;
        virtual basic_pheno get_pheno2(std::size_t) const;  
        
    protected:
        Phenotype(const_pheno_type); // these ones should not be available to the outside world
        Phenotype(pheno_type);       // 
        pheno_type pheno_ptr;
};

class PhenoBase    // This is the base class, it is just and interface (virtual pure)
{
	public:
		virtual ~PhenoBase() { }
        
        virtual pheno_type clone() const       = 0;
        virtual void add(const PhenoBase &)    = 0;
        virtual void remove(const PhenoBase &) = 0;
        virtual void divide(unsigned int)      = 0;
        virtual void square()                  = 0;   
        
        virtual std::size_t dimensionality() const = 0;
        virtual basic_pheno operator[] (std::size_t) const  = 0;
        virtual basic_pheno& operator[](std::size_t)        = 0;
        virtual basic_pheno get_pheno(std::size_t) const;
        virtual basic_pheno get_pheno2(std::size_t) const;      
        
        virtual void outformat(std::ostream&, unsigned int, unsigned int, std::string) const = 0;
        virtual void outformat2(std::ostream&, unsigned int, unsigned int, std::string) const;
        
	//~ private:
		//~ friend class boost::serialization::access;
		//~ template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            //~ ar & pheno;
            //~ ar & unstabpheno;
        //~ }  
};

class PhenoScalar: public PhenoBase {
  // The phenotype is a single scalar
  
  // Phenotypes should be canonical classes
  public: 
    PhenoScalar();
    PhenoScalar(const PhenoScalar&);
    PhenoScalar(basic_pheno);
    PhenoScalar& operator=(const PhenoScalar&);
    ~PhenoScalar() { }
    
    pheno_type clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    basic_pheno operator[] (std::size_t) const;
    basic_pheno& operator[](std::size_t);
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;    
          
  protected:
    
    PhenoScalar  operator +  (const PhenoScalar&) const;
    PhenoScalar& operator += (const PhenoScalar&);
    PhenoScalar  operator -  (const PhenoScalar&) const;
    PhenoScalar& operator -= (const PhenoScalar&);    
    PhenoScalar  operator /  (unsigned int) const;
    PhenoScalar& operator /= (unsigned int);  
  
    basic_pheno pheno;
    
};

class PhenoVector: public PhenoBase {
  // The phenotype is a vector
  
  public: 
    PhenoVector(); 
    PhenoVector(size_t);
    PhenoVector(const PhenoVector&);
    PhenoVector(const basic_container&);
    PhenoVector& operator = (const PhenoVector& p);
    ~PhenoVector() { }    
    
    pheno_type clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    basic_pheno operator[] (std::size_t) const;
    basic_pheno& operator[](std::size_t);
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;      
        
  protected:
    PhenoVector  operator +  (const PhenoVector&) const;
    PhenoVector& operator += (const PhenoVector&);
    PhenoVector  operator -  (const PhenoVector&) const;
    PhenoVector& operator -= (const PhenoVector&);    
    PhenoVector  operator /  (unsigned int) const;
    PhenoVector& operator /= (unsigned int);  
  
    basic_container pheno;
};

class PhenoTranscriptome: public PhenoBase {
  // The phenotype is a transcriptome
  //   - gene expression levels
  //   - gene expression stability
  
  public: 
    PhenoTranscriptome();
    PhenoTranscriptome(size_t);
    PhenoTranscriptome(const PhenoTranscriptome&);
    PhenoTranscriptome(const basic_container&, const basic_container&);
    PhenoTranscriptome& operator= (const PhenoTranscriptome&);
    ~PhenoTranscriptome() { }
    
    pheno_type clone() const;     
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    basic_pheno operator[] (std::size_t) const;
    basic_pheno& operator[](std::size_t); 
    basic_pheno get_pheno2(std::size_t) const;  
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const; 
    void outformat2(std::ostream&, unsigned int, unsigned int, std::string) const;             
    
  protected: 
     
    PhenoTranscriptome  operator +  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator += (const PhenoTranscriptome&); 
    PhenoTranscriptome  operator -  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator -= (const PhenoTranscriptome&);     
    PhenoTranscriptome  operator /  (unsigned int) const;
    PhenoTranscriptome& operator /= (unsigned int); 
  
    basic_container expression;
    basic_container stability;
};

#endif // PHENOTYPE_H_INCLUDED
