// Copyright 2020      Arnaud Le Rouzic    <lerouzic@egce.cnrs-gif.fr>

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <iostream>
#include <string>
#include <map>

#ifndef IOTAR_H_INCLUDED
#define IOTAR_H_INCLUDED

std::map<std::string, std::istream*> ifstream_from_tar(std::string);
std::string my_basename(std::string, char seperator = '/');

# endif // IOTAR_H
