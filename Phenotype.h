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

typedef double pheno_type;
typedef std::vector<Phenotype> pheno_container;

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
        typedef std::vector<pheno_type> basic_container;

    public:
        Phenotype();
        Phenotype(const Phenotype &);
        Phenotype(const PhenoBase &);
        ~Phenotype();
        Phenotype& operator=(const Phenotype &);
        
        // Specialized constructors, call directly a subclass (good idea or not?)
        Phenotype(pheno_type);  // PhenoScalar
        Phenotype(size_t);       // PhenoVector
        Phenotype(const basic_container&); // PhenoVector
        Phenotype(const basic_container&, const basic_container&); // PhenoTranscriptome
        
        static Phenotype sum  (const pheno_container&);
        static Phenotype sumsq(const pheno_container&);
        static Phenotype mean (const pheno_container&);        
        static Phenotype var  (const pheno_container&);

    // Necessary for compatibility with old code, should probably be avoided   
        bool is_defined() const;
        std::size_t dimensionality() const;
        pheno_type operator[] (std::size_t) const;
        pheno_type &operator[] (std::size_t);
        virtual pheno_type get_pheno(std::size_t) const;
        virtual pheno_type get_pheno2(std::size_t) const;  
        
    protected:
        Phenotype(std::unique_ptr<const PhenoBase>); // these ones should not be available to the outside world
        Phenotype(std::unique_ptr<PhenoBase>);       // 
        std::unique_ptr<PhenoBase> pheno_ptr;
};

class PhenoBase    // This is the base class, it is just and interface (virtual pure)
{
	public:
		virtual ~PhenoBase() { }
        
        virtual std::unique_ptr<PhenoBase> clone() const = 0;
        virtual void add(const PhenoBase &)    = 0;
        virtual void remove(const PhenoBase &) = 0;
        virtual void divide(unsigned int)      = 0;
        virtual void square()                  = 0;
        
        virtual std::size_t dimensionality() const = 0;
        virtual pheno_type operator[] (std::size_t) const  = 0;
        virtual pheno_type& operator[](std::size_t)        = 0;
        virtual pheno_type get_pheno(std::size_t) const;
        virtual pheno_type get_pheno2(std::size_t) const;      
        
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
    PhenoScalar(pheno_type);
    PhenoScalar& operator=(const PhenoScalar&);
    ~PhenoScalar() { }
    
    std::unique_ptr<PhenoBase> clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t);
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;    
          
  protected:
    
    PhenoScalar  operator +  (const PhenoScalar&) const;
    PhenoScalar& operator += (const PhenoScalar&);
    PhenoScalar  operator -  (const PhenoScalar&) const;
    PhenoScalar& operator -= (const PhenoScalar&);    
    PhenoScalar  operator /  (unsigned int) const;
    PhenoScalar& operator /= (unsigned int);  
  
    pheno_type pheno;
};

class PhenoVector: public PhenoBase {
  // The phenotype is a vector
  
  public: 
    PhenoVector(); 
    PhenoVector(size_t);
    PhenoVector(const PhenoVector&);
    PhenoVector(const Phenotype::basic_container&);
    PhenoVector& operator = (const PhenoVector& p);
    ~PhenoVector() { }    
    
    std::unique_ptr<PhenoBase> clone() const;
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t);
    
    void outformat(std::ostream&, unsigned int, unsigned int, std::string) const;      
        
  protected:
    PhenoVector  operator +  (const PhenoVector&) const;
    PhenoVector& operator += (const PhenoVector&);
    PhenoVector  operator -  (const PhenoVector&) const;
    PhenoVector& operator -= (const PhenoVector&);    
    PhenoVector  operator /  (unsigned int) const;
    PhenoVector& operator /= (unsigned int);  
  
    Phenotype::basic_container pheno;
};

class PhenoTranscriptome: public PhenoBase {
  // The phenotype is a transcriptome
  //   - gene expression levels
  //   - gene expression stability
  
  public: 
    PhenoTranscriptome();
    PhenoTranscriptome(size_t);
    PhenoTranscriptome(const PhenoTranscriptome&);
    PhenoTranscriptome(const Phenotype::basic_container&, const Phenotype::basic_container&);
    PhenoTranscriptome& operator= (const PhenoTranscriptome&);
    ~PhenoTranscriptome() { }
    
    std::unique_ptr<PhenoBase> clone() const;     
    void add(const PhenoBase&);
    void remove(const PhenoBase&);
    void divide(unsigned int);
    void square();
    
    std::size_t dimensionality() const;
    pheno_type  operator[](std::size_t) const;
    pheno_type& operator[](std::size_t); 
    pheno_type  get_pheno2(std::size_t) const;  
    
    void outformat (std::ostream&, unsigned int, unsigned int, std::string) const; 
    void outformat2(std::ostream&, unsigned int, unsigned int, std::string) const;             
    
  protected: 
     
    PhenoTranscriptome  operator +  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator += (const PhenoTranscriptome&); 
    PhenoTranscriptome  operator -  (const PhenoTranscriptome&) const;
    PhenoTranscriptome& operator -= (const PhenoTranscriptome&);     
    PhenoTranscriptome  operator /  (unsigned int) const;
    PhenoTranscriptome& operator /= (unsigned int); 
  
    Phenotype::basic_container expression;
    Phenotype::basic_container stability;
};

#endif // PHENOTYPE_H_INCLUDED
