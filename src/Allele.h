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

#include "types.h"
#include "Parameters.h"
#include "Mutation.h"

#include <vector>
#include <memory>

#ifdef SERIALIZATION_TEXT
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

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
        Allele() { }
        Allele(const std::vector<allele_type>);
	    Allele(const std::vector<allele_type>, const std::vector<std::string>);
	    Allele(const Allele &);
	    virtual ~Allele() { }
	
	    //operator overload
	    bool operator== (const Allele&) const;
	    bool operator!= (const Allele&) const;
	
	    //functions
	    unsigned int all_size() const;
		static std::vector<allele_type> combine_add(const Allele &, const Allele &);
		static std::vector<allele_type> combine_mean(const Allele &, const Allele &);
        static std::vector<allele_type> combine_OR(const Allele &, const Allele &);
		
		virtual std::shared_ptr<Allele> make_mutant_at_site(size_t site, const MutationModel&) const;
		virtual std::shared_ptr<Allele> make_mutant_random_site(const MutationModel&) const;
		virtual std::shared_ptr<Allele> make_mutant_all_sites(const MutationModel&) const;
        std::vector<allele_type> get_raw() const;
	
	protected :
	    std::vector<allele_type> allele;
	    std::vector<std::string> type_allele; // should be const, but const would make it difficult to serialize. 
        
	private:
        #ifdef SERIALIZATION_TEXT
		friend class boost::serialization::access;
		template<class Archive> void serialize(Archive & ar, const unsigned int version) {
            ar & allele;
            ar & type_allele;
        }
        #endif
};

#endif // ALLELE_H_INCLUDED
