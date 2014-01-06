#include "Architecture.h"
#include "ArchiAdditive.h"
#include "ArchiMultilinear.h"
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;



// constructors/destructor

Architecture::Architecture()
{
    assert(false); // The default constructor should never be used.
}


Architecture::Architecture(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


Architecture::Architecture(const ParameterSet& param)
    : gmap (param)
    , nloc (param.getpar(GENET_NBLOC) -> GetInt())
    , sall (param.getpar(GENET_ALLSIZE) -> GetInt())
    , mutrate (vector<double> (0))
    , mutsd (vector<double> (0))
{
    for (int i = 0; i < nloc; i++)
    {
        mutrate.push_back(param.getpar(GENET_MUTRATES)->GetDouble(i));
        mutsd.push_back(param.getpar(GENET_MUTSD)->GetDouble(i));
    }

}


// instance and initialization

Architecture* Architecture::instance = NULL;


void Architecture::initialize(const ParameterSet& param)
{
    if (Architecture::instance != NULL)
    {
        delete Architecture::instance;
        Architecture::instance = NULL;
    }
    //Architecture::instance = new Architecture(param);

    string type_archi = param.getpar(TYPE_ARCHI)->GetString();
    if (type_archi==AR_add)
    {
        Architecture::instance = new ArchiAdditive(param);
    }
    else if (type_archi==AR_mult)
    {
        Architecture::instance = new ArchiMultilinear(param);
    }
    else
    {
        assert("Wrong architecture type");
    }
}


Architecture* Architecture::Get()
{
    assert(Architecture::instance != NULL);
    return(Architecture::instance);
}


Architecture* Architecture::Get(const ParameterSet* param)
{
    if (param == NULL)
    {
        return(Architecture::Get());
    }
    else
    {
        return(Architecture::Get(*param));
    }
}


Architecture* Architecture::Get(const ParameterSet& param)
{
    if (Architecture::instance == NULL)
    {
        Architecture::instance = new Architecture(param);
    }
    return(Architecture::instance);
}


// operator overload
/*
ostream& operator << (ostream& out, const Architecture& archi)
{
    out << "=== Genetic map ===" << endl;
    out << archi.gmap;
    out << endl;

    out << "=== Mutation rates ===" << endl;
    for (int i = 0; i < archi.nb_loc(); i++)
    {
        out << "Loc" << i+1 << "\t" << archi.mutation_rate(i) << endl;
    }
    out << endl;
    return(out);
}
*/

// functions

int Architecture::nb_loc() const
{
    return nloc;
}


int Architecture::all_size() const
{
    return sall;
}


double Architecture::mutation_rate(int locus) const
{
    assert (locus >= 0);
    assert (locus < nloc);
    return(mutrate[locus]);
}


double Architecture::mutation_sd(int locus) const
{
    assert(locus >= 0);
    assert(locus < nloc);
    return(mutsd[locus]);
}


double Architecture::recombination_rate(int locus) const
{
    assert(locus >=0);
    assert(locus < nloc-1);
    return(gmap.recombination_rate(locus));
}


// for inheritance

double Architecture::phenotypic_value(const Genotype&) const
{
    assert(false);
}

