#ifndef MAIN_H_INCLUDED
#define MAIN_H_INCLUDED


#include "Parameters.h"
#include "OutputFormat.h"
#include "Architecture.h"
#include "ArchiAdditive.h"
#include "ArchiMultilinear.h"
#include "GeneticMap.h"
#include "Allele.h"
#include "Haplotype.h"
#include "Genotype.h"
#include "Individual.h"
#include "Population.h"
#include "Environment.h"
#include "Fitness.h"

#include <iostream>
#include <string>


const std::string SIMUL_GENER = "SIMUL_GENER";
const std::string SIMUL_REPET = "SIMUL_REPET";
const std::string SIMUL_OUTPUT = "SIMUL_OUTPUT";
const std::string SIMUL_EXTRA = "SIMUL_EXTRA";
const std::string GENET_NBLOC = "GENET_NBLOC";
const std::string GENET_ALLSIZE = "GENET_ALLSIZE";
const std::string GENET_MUTRATES = "GENET_MUTRATES";
const std::string GENET_MUTSD    = "GENET_MUTSD";
const std::string GENET_RECRATES = "GENET_RECRATES";
const std::string GENET_EPSILON2 = "GENET_EPSILON2";
const std::string GENET_EPSILON3 = "GENET_EPSILON3";
const std::string INIT_PSIZE = "INIT_PSIZE";
const std::string INIT_ALLELES = "INIT_ALLELES";
const std::string INIT_ALL_FRQ = "INIT_ALL_FRQ";
const std::string ENVIRO_SD = "ENVIRO_SD";
const std::string FITNESS_TYPE = "FITNESS_TYPE";
const std::string FITNESS_STRENGTH = "FITNESS_STRENGTH";
const std::string FITNESS_OPTIMUM = "FITNESS_OPTIMUM";
const std::string FITNESS_FLUCT = "FITNESS_FLUCT";
const std::string FITNESS_STRENGTH2 = "FITNESS_STRENGTH2";
const std::string FITNESS_OPTIMUM2 = "FITNESS_OPTIMUM2";
const std::string FITNESS_PERIOD = "FITNESS_PERIOD";
const std::string TYPE_ARCHI = "TYPE_ARCHI";



#endif // MAIN_H_INCLUDED
