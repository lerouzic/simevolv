#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED


#include "Random.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cassert>


class Parameter
{
public:
    //constructors/destructor
    virtual ~Parameter() {}

    //input/output
    virtual void read(std::istream &) {}
    virtual void write(std::ostream &) const {}

    //functions
    virtual long int GetInt() const {assert(false && "function unavailable");}
    virtual double GetDouble() const {assert(false && "function unavailable");}
    virtual double GetDouble(int) const {assert(false && "function unavailable");}
    virtual std::string GetString() const {assert(false && "function unavailable");}
    virtual bool is_nil() const {assert(false && "function unavailable");}
};


class Parameter_int: public Parameter
{
public:
    //constructors/destructor
    Parameter_int(long int minimum=0, long int maximum=999);
    ~Parameter_int() {}

    //input/output
    void read(std::istream&);
    void write(std::ostream&) const;

    //functions
    long int Get() const;
    long int GetInt() const {return(Get());}
    double GetDouble() const {return(double(Get()));}
    void Set(long int);
    bool is_nil() const {return(value == 0);}

protected:
    long int value;
    long int min;
    long int max;
    bool initialized;
};


class Parameter_double: public Parameter
{
public:
    //constructors/destructor
    Parameter_double(double minimum=-999.9, double maximum=999.9);
    ~Parameter_double() {}

    //input/output
    void read(std::istream&);
    void write(std::ostream&) const;

    //functions
    double Get() const;
    long int GetInt() const {return(int(Get()));}
    double GetDouble() const {return(Get());}
    void Set(double);
    bool is_nil() const {return(value == 0.0);}

protected:
    double value;
    double min;
    double max;
    bool initialized;
};


class Parameter_vector_double: public Parameter
{
public:
    //constructors/destructor
    Parameter_vector_double(double minimum=-999.9, double maximum=999.9);
    ~Parameter_vector_double() {}

    //input/output
    void read(std::istream&);
    void write(std::ostream&) const;

    //functions
    std::vector<double> Get() const;
    double Get_element(int) const;
    double GetDouble(int el) const {return(Get_element(el));}
    void Set(const std::vector<double>&);
    void Add(double);

protected:
    std::vector<double> value;
    double min;
    double max;
    bool initialized;
};


class Parameter_gaussian: public Parameter
{
public:
    //constructors/destructor
    Parameter_gaussian(double minimum_mean=-999.9, double maximum_mean=999.9, double maximum_sd=999.9);
    ~Parameter_gaussian() {}

    //input/output
    void read(std::istream&);
    void write(std::ostream&) const;

    //functions
    void SetMean(double);
    void SetSd(double);
    double draw() const;
    double GetDouble() const {return(draw());}
    bool is_nil() const {return((mean.Get() == 0.0) && (sd.Get() == 0.0));}

protected:
    Parameter_double mean;
    Parameter_double sd;
};


class Parameter_string: public Parameter 
{
public:
    //constructors/destructor
    Parameter_string(const std::vector<std::string>);
    ~Parameter_string();
    
    //input/output
    void read(std::istream&);
    void write(std::ostream&) const;
    
    void Set(std::string);
    std::string GetString() const;
    
protected:
	std::vector<std::string> possible_values; 
	std::string value;
	bool initialized;
};




class ParameterSet
{
public:
    // constructors/destructor
    ParameterSet();
    ParameterSet(const std::string& file);
    ~ParameterSet();

    // initialization
    void initialize();

    //input/output
    void read(const std::string&);
    void write(std::ostream&) const;

    //function
    const Parameter * getpar(const std::string&) const;

protected:
    std::map<std::string, Parameter*> parameters;
};


#endif // PARAMETERS_H_INCLUDED
