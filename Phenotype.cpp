#include "Phenotype.h"

using namespace std;



// constructors and destructors

Phenotype::Phenotype()
{
	// Nothing yet.
}


Phenotype::Phenotype(double init) 
{
	// only one phenotypic dimension!
	vector<double> vv;
	vv.push_back(init);
	initialize(vv);
}


Phenotype::Phenotype(const vector<double> & init) 
{
	initialize(init);
}


Phenotype::Phenotype(const Phenotype & templ)
{
	copy(templ);
}


Phenotype::~Phenotype() 
{
	
}


// initialize

void Phenotype::initialize(const vector<double> & init) 
{
	pheno = init;
}


// operator overload 

Phenotype & Phenotype::operator = (const Phenotype & templ) 
{
	if (this == &templ)
        return (*this);
    copy(templ);
    return(*this);
}


void Phenotype::copy(const Phenotype & templ) 
{
	pheno = templ.pheno;
}


// functions

double Phenotype::operator[] (const unsigned int index) const 
{
	return(pheno[index]);
}


unsigned int Phenotype::dimensionality() const 
{
	return(pheno.size());
}


void Phenotype::write_debug (ostream& out) const
{
	for (unsigned int i = 0; i < dimensionality(); i++) {
		out << pheno[i];
		if (i < dimensionality()-2)
			out << "/";
	}
	out << endl;
}


void Phenotype::write_simple (ostream& out) const
{
	write_debug(out);
}
