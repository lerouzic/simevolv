#ifndef PARCONST_H_INCLUDED
#define PARCONST_H_INCLUDED

#include <boost/assign/list_of.hpp> // Makes vector initialization easier


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

// Fitness types
const std::string FT_nosel = "no_selection";
const std::string FT_linear = "linear";
const std::string FT_expo = "exponential";
const std::string FT_gauss = "gaussian";
const std::string FT_quad = "quadratic";
const std::string FT_truncup = "trunc_up";
const std::string FT_truncdown = "trunc_down";
const std::string FT_concave = "concave";
const std::string FT_convex = "convex";
const std::vector<std::string> FT_options = boost::assign::list_of
		(FT_nosel)(FT_linear)(FT_expo)(FT_gauss)
		(FT_quad)(FT_truncup)(FT_truncdown)(FT_concave)(FT_convex);
		
// Fluctuation types
const std::string FF_nofluct = "no_fluctuation";
const std::string FF_smooth = "smooth";
const std::string FF_pflips = "periodic_flips";
const std::string FF_sflips = "stochastic_flips";
const std::string FF_brown = "brownian";
const std::string FF_white = "white_noise";
const std::vector<std::string> FF_options = boost::assign::list_of
		(FF_nofluct)(FF_smooth)(FF_pflips)(FF_sflips)(FF_brown)(FF_white);
		
// Architecture types
const std::string AR_add = "additive";
const std::string AR_mult = "multilinear";
const std::vector<std::string> AR_options = boost::assign::list_of
		(AR_add)(AR_mult);

#endif // PARCONST_H_INCLUDED
