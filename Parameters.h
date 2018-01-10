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

#include "types.h"

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
	    virtual float_type GetDouble() const {assert(false && "function unavailable");}
	    virtual float_type GetDouble(int) const {assert(false && "function unavailable");}
	    virtual std::string GetString() const {assert(false && "function unavailable");}
	    virtual std::vector<float_type> GetVectorDouble() const {assert(false && "function unavailable");}
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
	    float_type GetDouble() const;
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
	    Parameter_double(float_type minimum=-999.9, float_type maximum=999.9);
	    ~Parameter_double() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    float_type Get() const;
	    long int GetInt() const;
	    float_type GetDouble() const;
	    void Set(float_type);
	    bool is_nil() const;
	
	protected:
	    float_type value;
	    float_type min;
	    float_type max;
};

class Parameter_vector_double: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_vector_double(float_type minimum=-999.9, float_type maximum=999.9);
	    ~Parameter_vector_double() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    std::vector<float_type> Get() const;
	    std::vector<float_type> GetVectorDouble() const;
	    float_type Get_element(int) const;
	    float_type GetDouble(int) const;
	    void Set(const std::vector<float_type>&);
	    void Add(float_type);
	
	protected:
	    std::vector<float_type> value;
	    float_type min;
	    float_type max;
};

class Parameter_gaussian: public Parameter
{
	public:
	    // constructors/destructor
	    Parameter_gaussian(float_type minimum_mean=-999.9, float_type maximum_mean=999.9, float_type maximum_sd=999.9);
	    ~Parameter_gaussian() {}
	
	    // input/output
	    void read(std::istream&);
	    void write(std::ostream&) const;
	
	    // functions
	    void SetMean(float_type);
	    void SetSd(float_type);
	    float_type draw() const;
	    float_type GetDouble() const;
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
