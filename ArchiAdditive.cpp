#include "ArchiAdditive.h"
#include "Parconst.h"
#include "main.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;



// constructors and destuctor

ArchiAdditive::ArchiAdditive()
{
    assert(false); // The default constructor should never be used.
}


ArchiAdditive::ArchiAdditive(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


ArchiAdditive::ArchiAdditive(const ParameterSet& param)
    : Architecture(param)
{
    for (int i = 0; i < nloc; i++)
    {
        mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i));
        mutsd.push_back(param.getpar(GENET_MUTSD)->GetDouble(i));
    }
}


// operator overload
/*
ostream& operator << (ostream& out, const ArchiAdditive& archi)
{
    out << "=== Type of model ===" << endl;
    out << endl << "Addititive model" << endl;
    out << endl;
    out << endl;

    out << "=== Genetic map ===" << endl;
    out << archi.gmap;
    out << endl;

    out << "=== Mutation rates ===" << endl;
    out << endl;
    for (int i = 0; i < archi.nb_loc(); i++)
    {
        out << "Loc" << i+1 << "\t" << archi.mutation_rate(i) << endl;
    }
    out << endl;
    return(out);
}
*/


// functions

double ArchiAdditive::phenotypic_value (const Genotype& genotype) const
{
    int nloc = nb_loc();
    int sall = all_size();
    vector<double> sumall_father(nloc);
    vector<double> sumall_mother(nloc);
    vector<double> sumloc(nloc);
    double phenotype=0.0;

    for (int loc = 0 ; loc < nloc ; loc++)
    {
        for (int all = 0 ; all < sall ; all++)
        {
            sumall_father[loc] += genotype.gam_father.haplotype[loc].allele[all];
            sumall_mother[loc] += genotype.gam_mother.haplotype[loc].allele[all];

            sumloc[loc] = sumall_father[loc] + sumall_mother[loc];
        }
    }

    for (int loc = 0; loc < nloc; loc++)
    {
        phenotype += sumloc[loc];
    }

    return(phenotype);
}

