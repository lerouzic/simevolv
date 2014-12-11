// Copyright 2007-2014 Arnaud Le Rouzic    <lerouzic@legs.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/



#ifndef OUTPUTFORMAT_H_INCLUDED
#define OUTPUTFORMAT_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <vector>



void outformat(std::ostream &, const double, unsigned int width=12, unsigned int precision=5, const std::string & sep="");
void outformat(std::ostream &, const int, unsigned int width=12, const std::string & sep="");

void outformat(std::ostream &, const std::vector<double> &, unsigned int width=12, unsigned int precision=5, const std::string & sep="");

void outformat(std::ostream &, const std::string &, unsigned int width=12, const std::string & sep="");
void outformat(std::ostream &, unsigned int, const std::string &, unsigned int width=12, const std::string & sep="");

#endif // OUTPUTFORMAT_H_INCLUDED
