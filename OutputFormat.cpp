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

