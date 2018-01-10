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

#include "types.h"
#include "PhenoBase.h"

#include <iostream>
#include <string>  
#include <vector>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unique_ptr.hpp>
#endif

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
        Phenotype(pheno_type);  // PhenoScalar
        Phenotype(size_t);       // PhenoVector
        Phenotype(const std::vector<pheno_type>&); // PhenoVector
        Phenotype(const std::vector<pheno_type>&, const std::vector<pheno_type>&); // PhenoTranscriptome
        
        static Phenotype sum  (const std::vector<Phenotype> &);
        static Phenotype sumsq(const std::vector<Phenotype> &);
        static Phenotype mean (const std::vector<Phenotype> &);        
        static Phenotype var  (const std::vector<Phenotype> &);

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
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            // Why the following is necessary remains a mystery (BOOST_CLASS_EXPORT does not seem to do its job)
            ar.template register_type<PhenoScalar>();
            ar.template register_type<PhenoVector>();
            ar.template register_type<PhenoTranscriptome>();
            
            ar & pheno_ptr;
        }
        #endif
};

void outformat(std::ostream &, const Phenotype &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");
void outformat2(std::ostream &, const Phenotype &,
                      unsigned int width=13, unsigned int precision=5, std::string sep="");

#endif // PHENOTYPE_H_INCLUDED
