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



#include "OutputFormat.h"

using namespace std;



// constructors and destructor

OutputFormat::OutputFormat()
    : debugStream(NULL)
    , xmlStream(NULL)
    , simpleStream(NULL)
    , summaryStream(NULL)
{
    null = new nullstream();
    debugStream = null;
    xmlStream   = null;
    simpleStream= null;
    summaryStream=null;
}


OutputFormat::~OutputFormat()
{
    delete null;
}


// instance and initialization

OutputFormat * OutputFormat::instance = NULL;


bool OutputFormat::isInitialized()
{
    return(OutputFormat::instance != NULL);
}


// functions (output)

void OutputFormat::SetDebug(ostream & debug)
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    OutputFormat::instance->debugStream = &debug;
}


void OutputFormat::SetXml(ostream & xml)
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    OutputFormat::instance->xmlStream = &xml;
}


void OutputFormat::SetSimple(ostream & simple)
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    OutputFormat::instance->simpleStream = &simple;
}


void OutputFormat::SetSummary(ostream & summary)
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    OutputFormat::instance->summaryStream = &summary;
}


ostream & OutputFormat::GetDebug()
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    return(*OutputFormat::instance->debugStream);
}


ostream & OutputFormat::GetXml()
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    return(*OutputFormat::instance->xmlStream);
}


ostream & OutputFormat::GetSimple()
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    return(*OutputFormat::instance->simpleStream);
}


ostream & OutputFormat::GetSummary()
{
    if (OutputFormat::instance == NULL)
        OutputFormat::instance = new OutputFormat();

    return(*OutputFormat::instance->summaryStream);
}

