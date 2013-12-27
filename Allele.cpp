#include "Allele.h"
#include "Architecture.h"
#include "Random.h"

using namespace std;



// constructors/destructor

Allele::Allele(int sall)
{
    sall = Allele::all_size();

    for(int i = 0; i < sall; i++)
    {
        allele.push_back(Random::randnum());
    }
}


// operator overload

int Allele::operator==(const Allele& other) const
{
    return((*this).allele == other.allele);
}


int Allele::operator!=(const Allele& other) const
{
    return(!(*this == other));
}


// functions

int Allele::all_size() const
{
    Architecture * archi = Architecture::Get();
    int sall = archi -> all_size();

    return sall;
}


void Allele::print() const
{
    //cout << "Allele :" << endl;
    for(unsigned int i = 0; i < allele.size(); i++)
    {
        cout << allele[i] << " ";
    }
    cout << endl;
}


void Allele::make_mutation(int loc)
{
    // A mutation affects randomly one of the "sites" of the allele

    Architecture * archi = Architecture::Get();

    int mutated_site = floor(all_size()*Random::randnum());
    double modifier = archi->mutation_sd(loc) * Random::randgauss();
    allele[mutated_site] += modifier;
}

