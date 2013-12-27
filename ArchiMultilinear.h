#ifndef ARCHIMULTILINEAR_H_INCLUDED
#define ARCHIMULTILINEAR_H_INCLUDED

#include "Architecture.h"



class ArchiMultilinear : public Architecture
{
public :
    //constructors/destructor
    ArchiMultilinear();
    ArchiMultilinear(const Architecture&);
    ArchiMultilinear(const ParameterSet&);
    ~ArchiMultilinear() {}

    // operator overload
    friend std::ostream& operator << (std::ostream&, const Architecture&);

    //functions
    double get_epsilon2(int, int) const;
    double get_epsilon3(int, int, int) const;
    void set_epsilon2(int, int, double);
    void set_epsilon3(int, int, int, double);
    std::string print_epsilon2() const;
    std::string print_epsilon3() const;

    bool is_epistasis() const {return((is_epistasis2()) || (is_epistasis3()));}
    bool is_epistasis2() const {return(flag_epistasis2);}
    bool is_epistasis3() const {return(flag_epistasis3);}

    double phenotypic_value(const Genotype&) const;

protected :
    std::vector<std::vector<double> > epsilon2;
    std::vector<std::vector<std::vector<double> > > epsilon3;
    bool flag_epistasis2;
    bool flag_epistasis3;

};


#endif // ARCHIMULTILINEAR_H_INCLUDED
