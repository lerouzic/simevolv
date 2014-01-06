#ifndef ALLELE_H_INCLUDED
#define ALLELE_H_INCLUDED

#include <vector>



class Allele
{
    friend class Haplotype;
    friend class Architecture;
    friend class ArchiAdditive;
    friend class ArchiMultilinear;
	
	public :
	    //constructors/destructor
	    Allele(int nall);
	
	    //operator overload
	    int operator== (const Allele&) const;
	    int operator!= (const Allele&) const;
	
	    //functions
	    int all_size() const;
	    void print() const;
	    void make_mutation(int);
	
	protected :
	    std::vector<double> allele;
};


#endif // ALLELE_H_INCLUDED
