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



#ifndef PARCONST_H_INCLUDED
#define PARCONST_H_INCLUDED

// Boost : makes vector initialization easier
#include <boost/assign/list_of.hpp> 
#include <vector>



// COMMON PARAMETERS
const std::string TYPE_ARCHI = "TYPE_ARCHI";				/* Architecture types */
const std::string SIMUL_GENER = "SIMUL_GENER"; 
const std::string SIMUL_OUTPUT = "SIMUL_OUTPUT";
const std::string GENET_NBLOC = "GENET_NBLOC";
const std::string GENET_MUTRATES = "GENET_MUTRATES";
const std::string GENET_MUTSD    = "GENET_MUTSD";
const std::string GENET_RECRATES = "GENET_RECRATES";
const std::string INIT_PSIZE = "INIT_PSIZE";
const std::string INIT_ALLELES = "INIT_ALLELES";
const std::string INIT_CLONAL = "INIT_CLONAL";	 			/* Initial clonal status */
const std::string ENVIRO_SD = "ENVIRO_SD";
const std::string FITNESS_TYPE = "FITNESS_TYPE";			/* Selection types */
const std::string FITNESS_STRENGTH = "FITNESS_STRENGTH";
const std::string FITNESS_OPTIMUM = "FITNESS_OPTIMUM";
const std::string FITNESS_FLUCT = "FITNESS_FLUCT";			/* Fluctuation types */
const std::string FITNESS_STRENGTH2 = "FITNESS_STRENGTH2";
const std::string FITNESS_OPTIMUM2 = "FITNESS_OPTIMUM2";
const std::string FITNESS_PERIOD = "FITNESS_PERIOD";
const std::string FITNESS_STAB = "FITNESS_STAB";            /* Type of stability selection */
const std::string FITNESS_STABSTR = "FITNESS_STABSTR";      /* Strenght of selection on stability */
const std::string OUT_CANAL_TESTS = "OUT_CANAL_TESTS";
const std::string OUT_HERIT_TESTS = "OUT_HERIT_TESTS";
const std::string OUT_DIREPI_TESTS= "OUT_DIREPI_TESTS";

// MULTILINEAR ARCHITECTURE
const std::string GENET_EPSILON2 = "GENET_EPSILON2";
const std::string GENET_EPSILON3 = "GENET_EPSILON3";

// REGULATORY ARCHITECTURE
const std::string INIT_CONNECT = "INIT_CONNECT";
const std::string TYPE_SO = "TYPE_SO";	
const std::string INIT_BASAL = "INIT_BASAL";				/* Initial vector (So) types */
const std::string DEV_TIMESTEPS = "DEV_TIMESTEPS";
const std::string DEV_CALCSTEPS = "DEV_CALCSTEPS";




// Initial clonal status
const std::string CL_clonal = "clonal";
const std::string CL_variable = "variable";
const std::vector<std::string> CL_options = boost::assign::list_of (CL_clonal)(CL_variable);

// Selection types
const std::string FT_nosel = "no_selection";
const std::string FT_linear = "linear";
const std::string FT_expo = "exponential";
const std::string FT_gauss = "gaussian";
const std::string FT_quad = "quadratic";
const std::string FT_truncup = "trunc_up";
const std::string FT_truncdown = "trunc_down";
const std::string FT_concave = "concave";
const std::string FT_convex = "convex";
const std::vector<std::string> FT_options = boost::assign::list_of (FT_nosel)(FT_linear)(FT_expo)(FT_gauss)(FT_quad)(FT_truncup)(FT_truncdown)(FT_concave)(FT_convex);

// Fluctuation types
const std::string FF_nofluct = "no_fluctuation";
const std::string FF_smooth = "smooth";
const std::string FF_pflips = "periodic_flips";
const std::string FF_sflips = "stochastic_flips";
const std::string FF_brown = "brownian";
const std::string FF_white = "white_noise";
const std::vector<std::string> FF_options = boost::assign::list_of (FF_nofluct)(FF_smooth)(FF_pflips)(FF_sflips)(FF_brown)(FF_white);

// Strenght of selection on stability
const std::string FS_nostab = "no_stabsel";
const std::string FS_expo = "exponential_stab";
const std::vector<std::string> FS_options = boost::assign::list_of (FS_nostab)(FS_expo);

// Architecture types
const std::string AR_add = "additive";
const std::string AR_mult = "multilinear";
const std::string AR_wagner = "wagner";
const std::string AR_masel = "masel";
const std::string AR_siegal = "siegal";
const std::string AR_m2 = "m2";
const std::vector<std::string> AR_options = boost::assign::list_of (AR_add)(AR_mult)(AR_wagner)(AR_masel)(AR_siegal)(AR_m2);

// Initial vector type
const std::string SO_min = "minimum";
const std::string SO_max = "maximum";
const std::string SO_med = "median";
const std::string SO_randbin = "random_binary";
const std::string SO_rand = "random";
const std::string SO_basal = "basal";
const std::vector<std::string> SO_options = boost::assign::list_of (SO_min)(SO_max)(SO_med)(SO_randbin)(SO_rand)(SO_basal);


#endif // PARCONST_H_INCLUDED
