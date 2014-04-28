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



#ifndef OUTPUTFORMAT_H_INCLUDED
#define OUTPUTFORMAT_H_INCLUDED

#include <iostream>



struct nullstream:
std::ostream
{
    nullstream(): std::ios(0), std::ostream(0) {}
};


class OutputFormat
	{
	public:
	    // constructors/deconstructor
	    OutputFormat();
	    ~OutputFormat();
	
	    // instance/initialization
	    static OutputFormat * instance;
	    static bool isInitialized();
	
	    // output
	    static void SetDebug (std::ostream&);
	    static void SetXml   (std::ostream&);
	    static void SetSimple(std::ostream&);
	    static void SetSummary(std::ostream&);
	    static std::ostream& GetDebug();
	    static std::ostream& GetXml();
	    static std::ostream& GetSimple();
	    static std::ostream& GetSummary();
	
	protected:
	    std::ostream * null;
	    std::ostream * debugStream;
	    std::ostream * xmlStream;
	    std::ostream * simpleStream;
	    std::ostream * summaryStream;

};


#endif // OUTPUTFORMAT_H_INCLUDED
