#include "Parameters.h"
#include "Parconst.h"
#include "Random.h"
#include "main.h"

#include <algorithm> // for function std::find()

using namespace std;



// Parameter and daughter classes //

Parameter_int::Parameter_int(long int minimum, long int maximum)
    : min(minimum)
    , max(maximum)
    , initialized(false)
{}


long int Parameter_int::Get() const
{
    assert (initialized && "Parameter not initialized");
    return(value);
}


void Parameter_int::Set(long int v)
{
    initialized = true;
    assert (v >= min && v <= max && "Parameter out of range");
    value = v;
}


void Parameter_int::read(istream & i)
{
    long int val;
    assert(i >> val && "Unable to read parameter");
    Set(val);
}


void Parameter_int::write(ostream & out) const
{
    if (initialized)
    {
        out << value;
    }
    else
    {
        out << "Integer (between " << min << " and " << max << ")";
    }
}




Parameter_double::Parameter_double(double minimum, double maximum)
    : min(minimum)
    , max(maximum)
    , initialized(false)
{}


double Parameter_double::Get() const
{
    assert (initialized && "Parameter not initialized");
    return(value);
}


void Parameter_double::Set(double v)
{
    initialized = true;
    assert (v >= min && v <= max && "Parameter out of range");
    value = v;
}


void Parameter_double::read(istream & i)
{
    double val;
    assert (i >> val && "Unable to read parameter");
    Set(val);
}


void Parameter_double::write(ostream & out) const
{
    if (initialized)
    {
        out << value;
    }
    else
    {
        out << "Real (between " << min << " and " << max << ")";
    }
}




Parameter_vector_double::Parameter_vector_double(double minimum, double maximum)
    : min(minimum)
    , max(maximum)
    , initialized(false)
{
}


std::vector<double> Parameter_vector_double::Get() const
{
    assert (initialized && "Parameter not initialized");
    return(value);
}


double Parameter_vector_double::Get_element(int elem) const
{
    assert (initialized && "Parameter not initialized");
    if (elem < int(value.size()))
        return(value[elem]);
    else
        return(value[0]); // A warning could be emitted there
}


void Parameter_vector_double::Set(const vector<double>& v)
{
    for (vector<double>::const_iterator it = v.begin();
            it != v.end(); it++)
    {
        assert(*it >= min && *it <= max && "Parameter out of range");
    }
    initialized=true;
    value = v;
}

void Parameter_vector_double::Add(double e)
{
    assert (e >= min && e <= max && "Parameter out of range");
    initialized=true;
    value.push_back(e);
}


void Parameter_vector_double::read(istream & i)
{
    double val;
    while (i >> val)
    {
        Add(val);
    }
}


void Parameter_vector_double::write(ostream & out) const
{
    if (initialized)
    {
        for (vector<double>::const_iterator it = value.begin();
                it != value.end(); it++)
        {
            out << *it << "\t";
        }
    }
    else
    {
        out << "Real1 (between " << min << " and " << max << ")\t" << "Real2\tReal3\t...";
    }
}




Parameter_gaussian::Parameter_gaussian(double minimum_mean, double maximum_mean, double maximum_sd)
    : mean(minimum_mean, maximum_mean)
    , sd(0.0, maximum_sd)
{
}


void Parameter_gaussian::SetMean(double m)
{
    mean.Set(m);
}


void Parameter_gaussian::SetSd(double s)
{
    sd.Set(s);
}


double Parameter_gaussian::draw() const
{
    return(mean.Get() + sd.Get()*Random::randgauss());
}


void Parameter_gaussian::read(istream & i)
{
    double val;
    assert(i >> val &&  "Unable to read parameter");
    SetMean(val);
    assert(i >> val &&  "Unable to read parameter");
    SetSd(val);
}


void Parameter_gaussian::write(ostream & out) const
{
    mean.write(out);
    out << "\t";
    sd.write(out);
    out << "\t";
}


// Parameter_string

Parameter_string::Parameter_string(const vector<string> posval)
	:possible_values(posval)
	,initialized(false)
{
}


Parameter_string::~Parameter_string() 
{
}


void Parameter_string::read(istream & i)
{
    string val;
    assert (i >> val && "Unable to read parameter");
    Set(val);
}


void Parameter_string::write(ostream & out) const
{
    if (initialized)
    {
        out << value;
    }
    else
    {
        out << "String among: ";
        for (unsigned int i = 0; i < possible_values.size(); i++) {
			out << possible_values[i] << ",";
		}
		out << endl;
    }
}


void Parameter_string::Set(string v) {
	if (find(possible_values.begin(), possible_values.end(), v) != possible_values.end()) {
		value = v;
		initialized=true;
	} else {
		cerr << "Option " << v << " not recognized." << endl;
		cerr << "Accetable options are :" << endl;
		for (unsigned int i = 0; i < possible_values.size(); i++) {
			cerr << "\t" << possible_values[i] << endl;
		}
		assert("Stopping."); // This is bad error handling!
	}
}


string Parameter_string::GetString() const {
    assert (initialized && "Parameter not initialized");
	return(value);
}




// ParameterSet //

ParameterSet::ParameterSet()
{
    initialize();
}


ParameterSet::ParameterSet(const string & file)
{
    initialize();
    read(file);
}


ParameterSet::~ParameterSet()
{
    for (map<string, Parameter*>::iterator it = parameters.begin();
            it != parameters.end(); it++)
    {
        delete (it->second);
    }
}


void ParameterSet::initialize()
{
    // Simulation parameters
    parameters[SIMUL_GENER] = new Parameter_int(1, 1000*1000);
    parameters[SIMUL_REPET] = new Parameter_int(1, 10*1000);
    parameters[SIMUL_OUTPUT] = new Parameter_int(1, 10*1000);
    parameters[SIMUL_EXTRA] = new Parameter_int(0, 10*1000);

    // General genetic parameters
    parameters[GENET_NBLOC] = new Parameter_int(1, 100);
    parameters[GENET_ALLSIZE] = new Parameter_int(1, 100);
    parameters[GENET_MUTRATES] = new Parameter_vector_double(0.0, 1.0);
    parameters[GENET_MUTSD] = new Parameter_vector_double(0.0, 999.9);
    parameters[GENET_RECRATES] = new Parameter_vector_double(0.0, 0.5);
    parameters[GENET_EPSILON2] = new Parameter_gaussian(-999.9, 999.9, 999.9);
    parameters[GENET_EPSILON3] = new Parameter_gaussian(-999.9, 999.9, 999.9);

    // Initial population parameters
    parameters[INIT_PSIZE] = new Parameter_int(1, 1000*1000);
    parameters[INIT_ALLELES] = new Parameter_gaussian(-999.9,999.9,999.9);

    // Environmental parameters
    parameters[ENVIRO_SD] = new Parameter_double(0.0, 999.9);

    // Fitness parameters
    
    parameters[FITNESS_TYPE] = new Parameter_string(FT_options);
    parameters[FITNESS_STRENGTH] = new Parameter_double(-1000.*1000., 1000.*1000.);
    parameters[FITNESS_OPTIMUM] = new Parameter_double(-999.9, 999.9);
    parameters[FITNESS_FLUCT] = new Parameter_string(FF_options);
    parameters[FITNESS_OPTIMUM2] = new Parameter_double(-999.9, 999.9);
    parameters[FITNESS_STRENGTH2] = new Parameter_double(-999.9, 999.9);
    parameters[FITNESS_PERIOD] = new Parameter_int(0,1000*100);

    // Architecture type
    parameters[TYPE_ARCHI] = new Parameter_string(AR_options);
}


void ParameterSet::write(ostream & out) const
{
    for(map<string, Parameter*>::const_iterator it = parameters.begin();
            it != parameters.end(); it++)
    {
        out << it->first << "\t";
        it->second->write(out);
        out << endl;
    }
}


void ParameterSet::read(const string & file)
{
    ifstream paramfile(file.c_str());
    assert(paramfile && "Problem opening the parameter file");

    string line;
    while(!getline(paramfile, line).eof())
    {
        istringstream l(line);
        string name;
        l >> name;
        assert(parameters.find(name)!=parameters.end() && "Parameter key unknown.");
        parameters[name]->read(l);
    }

}


const Parameter * ParameterSet::getpar(const string & tag) const
{
    const Parameter * ans = parameters.find(tag)->second;
    return(ans);
}
