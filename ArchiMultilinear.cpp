// Copyright 2004-2007 José Alvarez-Castro <jose.alvarez-castro@lcb.uu.se>
// Copyright 2007      Arnaud Le Rouzic    <a.p.s.lerouzic@bio.uio.no>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#include "ArchiMultilinear.h"
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

ArchiMultilinear::ArchiMultilinear()
{
    assert(false); // The default constructor should never be used.
}


ArchiMultilinear::ArchiMultilinear(const Architecture& archi)
{
    assert(false); // The copy constructor should never be used.
}


ArchiMultilinear::ArchiMultilinear(const ParameterSet& param)
    : Architecture(param)
    , epsilon2(vector<vector<double> >(0))
    , epsilon3(vector<vector<vector<double> > >(0))
{
	// mutrate and mutsd are already intialized in the constructor of the parent class
    for (int loc1 = 0; loc1 < nloc; loc1++)
    {
        if (nloc > 1)
        {
            for (int loc2 = loc1+1; loc2 < nloc; loc2++)
            {
                set_epsilon2(loc1, loc2, param.getpar(GENET_EPSILON2)->GetDouble());
                if (nloc > 2)
                {
                    for (int loc3 = loc2+1; loc3 < nloc; loc3++)
                    {
                        set_epsilon3(loc1, loc2, loc3, param.getpar(GENET_EPSILON3)->GetDouble());
                    }
                }
            }
        }
    }
    flag_epistasis2 = !param.getpar(GENET_EPSILON2)->is_nil();
    flag_epistasis3 = !param.getpar(GENET_EPSILON3)->is_nil();
}

//operator overload
/*
std::ostream& operator << (std::ostream& out, const ArchiMultilinear& archi)
{
    out << "=== Type of model ===" << endl;
    out << endl << "Multilinear model" << endl;
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

    out << "=== 2-order epistasis ===" << std::endl;
    out << archi.print_epsilon2();
    out << endl;
    out << "=== 3-order epistasis ===" << std::endl;
    out << archi.print_epsilon3();
    out << endl;
    return(out);
}*/


// functions

double ArchiMultilinear::get_epsilon2(int loc1, int loc2) const
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    //assert (loc1 != loc2);
    //assert (loc1 >= 0);
    //assert (loc2 < nloc());
    //assert(int(epsilon2.size()) >= loc1+1);
    //assert(int(epsilon2[loc1].size()) >= loc2-loc1);
    return(epsilon2[loc1][loc2-loc1-1]);
}


double ArchiMultilinear::get_epsilon3(int loc1, int loc2, int loc3) const
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    if (loc2 > loc3)
    {
        swap(loc2, loc3);
    }
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    //assert(loc1 != loc2);
    //assert(loc2 != loc3);
    //assert(loc1 >= 0);
    //assert(loc3 < nloc());
    //assert(int(epsilon3.size()) >= loc1+1);
    //assert(int(epsilon3[loc1].size()) >= loc2-loc1);
    //assert(int(epsilon3[loc1][loc2-loc1-1].size()) >= loc3-loc2);
    return(epsilon3[loc1][loc2-loc1-1][loc3-loc2-1]);
}


void ArchiMultilinear::set_epsilon2(int loc1, int loc2, double value)
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert (loc1 != loc2);
    assert (loc1 >= 0);
    assert (loc2 < nb_loc());

    while(int(epsilon2.size()) < (loc1+1))
    {
        vector<double> v(0);
        epsilon2.push_back(v);
    }

    while(int(epsilon2[loc1].size()) < (loc2-loc1))
    {
        epsilon2[loc1].push_back(0);
    }

    epsilon2[loc1][loc2-loc1-1] = value;
}


void ArchiMultilinear::set_epsilon3(int loc1, int loc2, int loc3, double value)
{
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    if (loc2 > loc3)
    {
        swap(loc2, loc3);
    }
    if (loc1 > loc2)
    {
        swap(loc1, loc2);
    }
    assert(loc1 != loc2);
    assert(loc2 != loc3);
    assert(loc1 >= 0);
    assert(loc3 < nb_loc());

    while(int(epsilon3.size()) < nb_loc()-2)
    {
        vector<vector<double> > v(0);
        epsilon3.push_back(v);
    }

    while(int(epsilon3[loc1].size()) < (nb_loc()-loc1-2))
    {
        vector<double> v(0);
        epsilon3[loc1].push_back(v);
    }

    while(int(epsilon3[loc1][loc2-loc1-1].size()) < (nb_loc()-loc2-1))
    {
        epsilon3[loc1][loc2-loc1-1].push_back(0);
    }

    epsilon3[loc1][loc2-loc1-1][loc3-loc2-1] = value;
}


string ArchiMultilinear::print_epsilon2() const
{
    ostringstream out;

    out << "\t";
    for (int loc2 = 1; loc2 < nb_loc(); loc2++)
        out << "Loc" << loc2 << "\t";
    out << endl;
    for (int loc1 = 0; loc1 < nb_loc()-1; loc1++)
    {
        out << "Loc" << loc1 << "\t";
        for (int i = 0; i < loc1; i++)
            out << "\t";
        for (int loc2 = loc1+1; loc2 < nb_loc(); loc2++)
        {
            out << setiosflags(ios::fixed) << setprecision(3) << get_epsilon2(loc1, loc2) << "\t";
        }
        out << endl;
    }
    return(out.str());
}


string ArchiMultilinear::print_epsilon3() const
{
    ostringstream out;

    for (int loc1 = 0; loc1 < nb_loc() -2; loc1++)
    {
        out << "* Locus " << loc1 << endl;
        out << "\t";
        for (int loc3 = loc1+2; loc3 < nb_loc(); loc3++)
            out << "Loc" << loc3 << "\t";
        out << endl;
        for (int loc2 = loc1+1; loc2 < nb_loc()-1; loc2++)
        {
            out << "Loc" << loc2 << "\t";
            for (int i = 0; i < loc2-loc1-1; i++)
                out << "\t";
            for (int loc3 = loc2+1; loc3 < nb_loc(); loc3++)
                out << setiosflags(ios::fixed) << setprecision(3) << get_epsilon3(loc1, loc2, loc3) << "\t";
            out << endl;
        }
        out << endl;
    }
    return(out.str());
}


Phenotype ArchiMultilinear::phenotypic_value (const Genotype& genotype) const
{
    int nloc = nb_loc();
    int sall = all_size();
    vector<double> sumall_father(nloc);
    vector<double> sumall_mother(nloc);
    vector<double> sumloc(nloc);
    double phenotype=0.0;

    for (int loc1 = 0 ; loc1 < nloc ; loc1++)
    {
        for (int all = 0 ; all < sall ; all++)
        {
            sumall_father[loc1] += genotype.gam_father.haplotype[loc1]->allele[all];
            sumall_mother[loc1] += genotype.gam_mother.haplotype[loc1]->allele[all];

            sumloc[loc1] = sumall_father[loc1] + sumall_mother[loc1];
        }
    }

    for (int loc1 = 0; loc1 < nloc; loc1++)
    {
        phenotype += sumloc[loc1];
        for (int loc2 = loc1+1; loc2 < nloc; loc2++)
        {
            if(is_epistasis2())
            {
                phenotype += get_epsilon2(loc1, loc2) * sumloc[loc1] * sumloc[loc2];
            }

            for (int loc3 = loc2+1; loc3 < nloc; loc3++)
            {
                if(is_epistasis3())
                {
                    phenotype += get_epsilon3(loc1, loc2, loc3) * sumloc[loc1] * sumloc[loc2] * sumloc[loc3];
                }
            }
        }
    }
    return(Phenotype(phenotype));
}

