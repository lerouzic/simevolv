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



#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cassert>
#include <map>


////////////////////////////////// PARAMETER ////////////////////////////////////////

class Parameter
{ // virtual
	public:
	    //constructors/destructor
	    Parameter();
	    virtual ~Parameter() {}
	
	    //input/output
	    virtual void read(std::istream &) {}
	    virtual void write(std::ostream &) const {}
	
	    //functions
	    virtual long int GetInt() const {assert(false && "function unavailable");}
	    virtual double GetDouble() const {assert(false && "function unavailable");}
	    virtual double GetDouble(int) const {assert(false && "function unavailable");}
	    virtual std::string GetString() const {assert(false && "function unavailable");}
	    virtual std::vector<double> GetVectorDouble() const {assert(false && "function unavailable");}
	    virtual bool is_nil() const {assert(false && "function unavailable");}
	    virtual bool is_initialized() const;
	    unsigned int get_count() const;
	    
	protected:
		mutable unsigned int count;
		bool initialized;
};

class Parameter_int: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_int(long int minimum=0, long int maximum=999);
	    ~Parameter_int() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    long int Get() const;
	    long int GetInt() const;
	    double GetDouble() const;
	    void Set(long int);
	    bool is_nil() const;
	
	protected:
	    long int value;
	    long int min;
	    long int max;
};

class Parameter_double: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_double(double minimum=-999.9, double maximum=999.9);
	    ~Parameter_double() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    double Get() const;
	    long int GetInt() const;
	    double GetDouble() const;
	    void Set(double);
	    bool is_nil() const;
	
	protected:
	    double value;
	    double min;
	    double max;
};

class Parameter_vector_double: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_vector_double(double minimum=-999.9, double maximum=999.9);
	    ~Parameter_vector_double() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    std::vector<double> Get() const;
	    std::vector<double> GetVectorDouble() const;
	    double Get_element(int) const;
	    double GetDouble(int) const;
	    void Set(const std::vector<double>&);
	    void Add(double);
	
	protected:
	    std::vector<double> value;
	    double min;
	    double max;
};

class Parameter_gaussian: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_gaussian(double minimum_mean=-999.9, double maximum_mean=999.9, double maximum_sd=999.9);
	    ~Parameter_gaussian() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    void SetMean(double);
	    void SetSd(double);
	    double draw() const;
	    double GetDouble() const;
	    bool is_nil() const;
	    bool is_initialized() const;
	
	protected:
	    Parameter_double mean;
	    Parameter_double sd;
};

class Parameter_string: public Parameter 
{
	public:
	    // constructors/destructor
	    Parameter_string();
	    Parameter_string(const std::vector<std::string>);
	    ~Parameter_string() {}
	    
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	    
	    // functions
	    void Set(std::string);
	    std::string GetString() const;
	    
	protected:
		std::vector<std::string> possible_values; 
		std::string value;
};



////////////////////////////////// PARAMETER SET ////////////////////////////////////////

class ParameterSet
{
	public:
	    // constructors/destructor
	    ParameterSet();
	    ParameterSet(const std::string& file);
	    ~ParameterSet();
	
	    // initialization
	    void initialize();
	
	    // input/output
	    void read(const std::string&);
	    void write(std::ostream&) const;
	
	    // function
	    const Parameter * getpar(const std::string&) const;
	    bool exists(const std::string&) const;
	    
	    // consistency checks
	    void warning_unused() const;
	    void warning_multicalls(unsigned int morethan = 1) const;
	
	protected:
	    std::map<std::string, Parameter*> parameters;
};

#endif // PARAMETERS_H_INCLUDED
