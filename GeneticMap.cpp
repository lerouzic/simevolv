#include "GeneticMap.h"
#include "Architecture.h"
#include "main.h"

using namespace std;


//constructors and destructor

GeneticMap::GeneticMap()
    : recrate(vector<double> (0))
{
}


GeneticMap::GeneticMap(const ParameterSet& param)
    : recrate(vector<double> (0))
{
    int nloc = param.getpar(GENET_NBLOC) -> GetInt();

    for (int i = 0; i < nloc-1; i++)
    {
        recrate.push_back(param.getpar(GENET_RECRATES)->GetDouble(i));
    }
}


// operator overload

std::ostream& operator << (std::ostream& out, const GeneticMap& gmap)
{
    int nloc = gmap.nb_loc();
    out << endl << "Number of loci : " << endl << nloc << endl;

    out << endl << "Recombination rates : " << endl;
    for (int i = 0; i < nloc-1; i++)
    {
        out << gmap.recrate[i] << endl;
    }
    out << endl;

    return(out);
}


//functions

int GeneticMap::nb_loc() const
{
    int nloc = (recrate.size()+1);

    return nloc;
}


double GeneticMap::recombination_rate(int loc1, int loc2) const
{
    int nloc = GeneticMap::nb_loc();

    if (loc2 == -1)
    {
        loc2 = loc1 + 1;
    }

    assert(loc1 >= 0);
    assert(loc2 < nloc);

    if (loc2 == loc1 + 1)
    {
        return(recrate[loc1]);
    }
    else
    {
        assert(1); // does not work yet (is it useful anyway?)
        double dist=0;
        for (int i = loc1; i < loc2-1; i++)
        {
            dist += recrate[loc1];
        }
        return(exp(-2*dist));
    }
}

