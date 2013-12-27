#ifndef ARCHIADDITIVE_H_INCLUDED
#define ARCHIADDITIVE_H_INCLUDED

#include "Architecture.h"



class ArchiAdditive : public Architecture
{
public :
    //constructors/destructor
    ArchiAdditive();
    ArchiAdditive(const Architecture&);
    ArchiAdditive(const ParameterSet&);
    ~ArchiAdditive() {}

    // operator overload
    friend std::ostream& operator << (std::ostream&, const Architecture&);

    //functions
    double phenotypic_value(const Genotype&) const;

protected :

};

#endif // ARCHIADDITIVE_H_INCLUDED
