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



#include "Parameters.h"

#include "Parconst.h"
#include "Random.h"

#include <algorithm> // for function std::find()

using namespace std;


///// Parameters

Parameter::Parameter()
	: count(0)
	, initialized(false)
{
}

bool Parameter::is_initialized() const 
{
	return(initialized);
}

unsigned int Parameter::get_count() const 
{
	return(count);
}


///// Daughter classes /////

Parameter_int::Parameter_int(long int minimum, long int maximum)
	: Parameter()
    , min(minimum)
    , max(maximum)
{
}

void Parameter_int::read(istream & i)
{
    long int val;
    assert(i >> val && "Unable to read parameter");
    Set(val);
}

void Parameter_int::write(ostream & out) const
{
    if (is_initialized())
    {
        out << value;
    }
    else
    {
        out << "Integer (between " << min << " and " << max << ")";
    }
}

void Parameter_int::erase() 
{
    initialized = false;
    value = 0; // there is no default, but better to set it to some expected value
}

long int Parameter_int::Get() const
{
    assert (is_initialized() && "Parameter not initialized");
    count++;
    return(value);
}

long int Parameter_int::GetInt() const 
{
	return(Get());
}

float_type Parameter_int::GetDouble() const 
{
	return(float_type(Get()));
}
	    
void Parameter_int::Set(long int v)
{
    initialized = true;
    assert (v >= min && v <= max && "Parameter out of range");
    value = v;
}

bool Parameter_int::is_nil() const 
{
	return(value == 0);
}


Parameter_double::Parameter_double(float_type minimum, float_type maximum)
    : Parameter()
    , min(minimum)
    , max(maximum)
{
}

void Parameter_double::read(istream & i)
{
    float_type val;
    assert (i >> val && "Unable to read parameter");
    Set(val);
}

void Parameter_double::write(ostream & out) const
{
    if (is_initialized())
    {
        out << value;
    }
    else
    {
        out << "Real (between " << min << " and " << max << ")";
    }
}

void Parameter_double::erase()
{
    initialized = false;
    value = 0.0; // no default, why not 0?
}

float_type Parameter_double::Get() const
{
    assert (is_initialized() && "Parameter not initialized");
    count++;
    return(value);
}

long int Parameter_double::GetInt() const 
{
	return(int(Get()));
}

float_type Parameter_double::GetDouble() const 
{
	return(Get());
}

void Parameter_double::Set(float_type v)
{
    initialized = true;
    assert (v >= min && v <= max && "Parameter out of range");
    value = v;
}

bool Parameter_double::is_nil() const 
{
	return(value == 0.0);
}


Parameter_vector_double::Parameter_vector_double(float_type minimum, float_type maximum)
    : Parameter()
    , min(minimum)
    , max(maximum)
{
}

void Parameter_vector_double::read(istream & i)
{
	value.clear();
    float_type val;
    while (i >> val)
    {
        Add(val);
    }
}

void Parameter_vector_double::write(ostream & out) const
{
    if (is_initialized())
    {
        for (vector<float_type>::const_iterator it = value.begin(); it != value.end(); it++)
        {
            out << *it << "\t";
        }
    }
    else
    {
        out << "Real1 (between " << min << " and " << max << ")\t" << "Real2\tReal3\t...";
    }
}

void Parameter_vector_double::erase()
{
    initialized = false;
    value.clear();
}

vector<float_type> Parameter_vector_double::Get() const
{
    assert (is_initialized() && "Parameter not initialized");
    count++;
    return(value);
}

vector<float_type> Parameter_vector_double::GetVectorDouble() const 
{
	return(Get());
}

float_type Parameter_vector_double::Get_element(int elem) const
{
    assert (is_initialized() && "Parameter not initialized");
    if (count == 0)
    {
		count++;  // otherwise accessing multiple elements counts several accesses
    }
    if (elem < int(value.size()))
    {
		return(value[elem]);
    }
    else
    {
        return(value[0]); // A warning could be emitted there
	}
}

float_type Parameter_vector_double:: GetDouble(int el) const 
{
	return(Get_element(el));
}

void Parameter_vector_double::Set(const vector<float_type>& v)
{
    for (vector<float_type>::const_iterator it = v.begin(); it != v.end(); it++)
    {
        assert(*it >= min && *it <= max && "Parameter out of range");
    }
    initialized=true;
    value = v;
}

void Parameter_vector_double::Add(float_type e)
{
    assert (e >= min && e <= max && "Parameter out of range");
    initialized=true;
    value.push_back(e);
}


Parameter_gaussian::Parameter_gaussian(float_type minimum_mean, float_type maximum_mean, float_type maximum_sd)
    : Parameter()
    , mean(minimum_mean, maximum_mean)
    , sd(0.0, maximum_sd)
{
}

void Parameter_gaussian::read(istream & i)
{
    float_type val;
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

void Parameter_gaussian::erase()
{
    initialized = false;
    mean = 0.0;
    sd   = 0.0; // no default values, why not 0?
}

void Parameter_gaussian::SetMean(float_type m)
{
    mean.Set(m);
}

void Parameter_gaussian::SetSd(float_type s)
{
    sd.Set(s);
}

float_type Parameter_gaussian::draw() const
{
	assert(is_initialized() && "Parameter not initialized");
	count ++;
    return(mean.Get() + sd.Get()*Random::randgauss());
}

float_type Parameter_gaussian::GetDouble() const 
{
	return(draw());
}

bool Parameter_gaussian::is_nil() const 
{
	return((mean.Get() == 0.0) && (sd.Get() == 0.0));
}

bool Parameter_gaussian::is_initialized() const 
{
	return(mean.is_initialized() && sd.is_initialized());
}
	


Parameter_string::Parameter_string()
	: Parameter()
	, possible_values()
	, value("NotInitialized")
{
}

Parameter_string::Parameter_string(const vector<string> posval)
	: Parameter()
	, possible_values(posval)
	, value("NotInitialized")
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
    if (is_initialized())
    {
        out << value;
    }
    else
    {
        out << "String among: ";
        for (unsigned int i = 0; i < possible_values.size(); i++)
        {
			out << possible_values[i] << ",";
		}
		out << endl;
    }
}

void Parameter_string::erase()
{
    initialized = false;
    value = "";
}

void Parameter_string::Set(string v) 
{
	if ((possible_values.empty()) || 
		find(possible_values.begin(), possible_values.end(), v) != possible_values.end()) 
	{
		value = v;
		initialized=true;
	} 
	else 
	{
		cerr << "Option " << v << " not recognized." << endl;
		cerr << "Acceptable options are :" << endl;
		for (unsigned int i = 0; i < possible_values.size(); i++) 
		{
			cerr << "\t" << possible_values[i] << endl;
		}
		exit(EXIT_FAILURE);
	}
}

string Parameter_string::GetString() const 
{
    assert (is_initialized() && "Parameter not initialized");
    count++;
	return(value);
}



///////// ParameterSet 

ParameterSet::ParameterSet()
{
    initialize();
}

ParameterSet::ParameterSet(const string & file)
{
    initialize();
    read(file);
}

ParameterSet::ParameterSet(istream * is)
{	
	initialize();
	read(is);
}

ParameterSet::~ParameterSet()
{
    for (auto it = parameters.begin();
            it != parameters.end(); it++)
    {
        delete (it->second);
    }
}

void ParameterSet::initialize()
{
    // Architecture type
    parameters[TYPE_ARCHI] = new Parameter_string(AR_options);
    
    // Simulation parameters
    parameters[SIMUL_GENER] = new Parameter_int(0, 1000*1000);
    parameters[SIMUL_MAXGEN] = new Parameter_int(0, 1000*1000);
    parameters[SIMUL_OUTPUT] = new Parameter_int(1, 10*1000);
   
	// General genetic parameters
    parameters[GENET_NBLOC] = new Parameter_int(1, 100);
    parameters[GENET_NBPHEN] = new Parameter_int(1, 100);
    parameters[GENET_PLOIDY] = new Parameter_int(1, 2);
    parameters[GENET_MUTTYPE] = new Parameter_string(MT_options);
    parameters[GENET_MUTMEM]  = new Parameter_string(MM_options);
    parameters[GENET_MUTRATES] = new Parameter_vector_double(0.0, 1.0);
    parameters[GENET_MUTSD] = new Parameter_vector_double(0.0, 999.9);
    parameters[GENET_RECRATES] = new Parameter_vector_double(0.0, 0.5);
    parameters[GENET_SELFING] = new Parameter_double(0.0, 1.0);
    parameters[GENET_CLONAL]  = new Parameter_double(0.0, 1.0);
    parameters[GENET_EPIGENET] = new Parameter_double(0.0, 1.0);
   
    // Initial population parameters
    parameters[INIT_PSIZE] = new Parameter_int(1, 1000*1000);
    parameters[INIT_ALLELES] = new Parameter_gaussian(-999.9,999.9,999.9);
    parameters[TYPE_ALLELES] = new Parameter_string(TA_options);
    parameters[INIT_CLONAL] = new Parameter_string(CL_options);

    // Environmental parameters
    parameters[ENVIRO_SDINIT] = new Parameter_double(0.0, 999.9);
	parameters[ENVIRO_SDDYNAM] = new Parameter_double(0.0, 999.9);
	parameters[ENVIRO_SDFINAL] = new Parameter_double(0.0, 999.9);
    parameters[ENVIRO_PLASTICITY] = new Parameter_vector_double(0.0, 1.0);

    // Fitness parameters
    parameters[FITNESS_TYPE] = new Parameter_string(FT_options);
    parameters[FITNESS_STRENGTH] = new Parameter_vector_double(-1000.*1000., 1000.*1000.);
    parameters[FITNESS_OPTIMUM] = new Parameter_vector_double(-999.9, 999.9);
    parameters[FITNESS_CORRELATION] = new Parameter_vector_double(-1.0, 1.0);
    
    // Output measurements
    parameters[OUT_GENO] = new Parameter_string(OG_options);
    parameters[OUT_UNSTAB] = new Parameter_string(OU_options);
    parameters[OUT_CANAL_TESTS] = new Parameter_int(0, 1000*1000);
    parameters[OUT_CANAL_MUTSD] = new Parameter_vector_double(0.0, 999.9);
    parameters[OUT_CANAL_SDINIT] = new Parameter_double(0.0, 999.9);
    parameters[OUT_CANAL_SDDYNAM] = new Parameter_double(0.0, 999.9);    
    parameters[OUT_HERIT_TESTS] = new Parameter_int(0, 1000*1000);
    parameters[OUT_DIREPI_TESTS] = new Parameter_int(0, 1000*1000);
    
    parameters[PHENO_SCALING] = new Parameter_string(ST_options);
        
    // Multilinear architecture
    parameters[GENET_EPSILON2e] = new Parameter_gaussian(-999.9, 999.9, 999.9);
    parameters[GENET_EPSILON2p] = new Parameter_gaussian(-999.9, 999.9, 999.9);
    //parameters[GENET_EPSILON3] = new Parameter_gaussian(-999.9, 999.9, 999.9);
    
    //Boolean architecture
    parameters[MATRIX_DENS] = new Parameter_double(0.0, 1.0);
    parameters[LOG_OPERATOR_DENS] = new Parameter_double(0.0, 1.0);
    parameters[PHEN_NBLOC] = new Parameter_int(1, 100);
    parameters[SCALE] = new Parameter_string(SC_options);
    
    // Regulatory architecture
    parameters[INIT_CONNECT] = new Parameter_double(0.0, 1.0);
    parameters[INIT_CONDIAG] = new Parameter_double(0.0, 1.0);
	parameters[TYPE_SO] = new Parameter_string(SO_options);
	parameters[INIT_BASAL] = new Parameter_double(0.0, 1.0);
	parameters[INIT_RECURRENCE] = new Parameter_double(0.0, 1.0);
	parameters[DEV_TIMESTEPS] = new Parameter_int(0, 100*1000);
	parameters[DEV_CALCSTEPS] = new Parameter_int(0, 100*1000);
	parameters[FITNESS_STAB] = new Parameter_string(FS_options);
    parameters[FITNESS_STABSTR] = new Parameter_vector_double(0.0, 1000.*1000.);
    
    // Input/Output
    parameters[FILE_NEXTPAR] = new Parameter_string();
}

void ParameterSet::write(ostream & out) const
{
    for(map<string, Parameter*>::const_iterator it = parameters.begin(); it != parameters.end(); it++)
    {
        out << it->first << "\t";
        it->second->write(out);
        out << endl;
    }
}

void ParameterSet::read(const string & file)
{
    ifstream paramstream(file.c_str());
    if (!paramstream) 
    {
		cerr << "Problem opening the parameter file " << file << endl;
		exit(EXIT_FAILURE);
	}
	read(&paramstream);
}

void ParameterSet::read(istream * paramstream) {
    string line;
    while(!getline(*paramstream, line).eof())
    {
        istringstream l(line);
        string name;
        l >> name;
        if (parameters.find(name)==parameters.end() ) 
        {
			cerr << "Parameter " << name << " unknown." << endl;
			exit(EXIT_FAILURE); 			
		}
        parameters[name]->read(l);
    }
}

void ParameterSet::erase(const string & tag)
{
	if (!exists(tag)) 
	{
		cerr << "Parameter " << tag << " is missing." << endl;
		exit(EXIT_FAILURE); 
	}
    auto pp = parameters.find(tag);
    pp->second->erase();  
}

const Parameter * ParameterSet::getpar(const string & tag) const
{
	if (!exists(tag)) 
	{
		cerr << "Parameter " << tag << " is missing." << endl;
		exit(EXIT_FAILURE); 
	}
    const Parameter * ans = parameters.find(tag)->second;
    
    return(ans);
}

bool ParameterSet::exists(const string & tag) const
{
	return((parameters.count(tag) != 0) && (parameters.find(tag)->second->is_initialized()));
}

void ParameterSet::warning_unused() const 
{
    for (auto it = parameters.begin(); it != parameters.end(); it++)
    {
		if (it->second->is_initialized() && it->second->get_count() == 0) 
		{
			cerr << "Warning: Parameter " << it->first << " initialized but unused." << endl;
		}
	}	
}

void ParameterSet::warning_multicalls(unsigned int morethan /* = 1 */) const
{
    for (auto it = parameters.begin(); it != parameters.end(); it++)
    {
		if (it->second->get_count() > morethan) 
		{
			cerr << "Warning: Parameter " << it->first << " called " << it->second->get_count() << " times." << endl;
		}		
	}	
}
