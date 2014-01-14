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



#include "Allele.h"
#include "Architecture.h"
#include "Random.h"

#include <iostream>
#include <string>
#include <cmath>

using namespace std;



// constructors and destructor

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

