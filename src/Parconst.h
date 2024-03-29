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

#include <boost/assign/list_of.hpp> // Boost : makes vector initialization easier
#include <vector>

const double MIN_LOG_VAR = -20.0;

// COMMON PARAMETERS
const std::string TYPE_ARCHI = "TYPE_ARCHI";				/* Architecture types */
const std::string SIMUL_GENER = "SIMUL_GENER"; 
const std::string SIMUL_OUTPUT = "SIMUL_OUTPUT";
const std::string SIMUL_MAXGEN = "SIMUL_MAXGEN";
const std::string GENET_NBLOC = "GENET_NBLOC";
const std::string GENET_NBPHEN = "GENET_NBPHEN";
const std::string GENET_MUTTYPE = "GENET_MUTTYPE"; 		/* Type of mutations: "individual" or "locus" */
const std::string GENET_MUTMEM  = "GENET_MUTMEM";		/* Memory of mutational effects: "cumul" or "stationary" */
const std::string GENET_MUTRATES = "GENET_MUTRATES";
const std::string GENET_PLOIDY = "GENET_PLOIDY";
const std::string GENET_MUTSD    = "GENET_MUTSD";
const std::string GENET_RECRATES = "GENET_RECRATES";
const std::string GENET_SELFING = "GENET_SELFING";
const std::string GENET_CLONAL  = "GENET_CLONAL";
const std::string GENET_EPIGENET = "GENET_EPIGENET";
const std::string INIT_PSIZE = "INIT_PSIZE";
const std::string INIT_ALLELES = "INIT_ALLELES";            /* Mean and variance of initial alleles */
const std::string INIT_ALLELES_FULL = "INIT_ALLELES_FULL";       /* Exhaustive vector of initial allelic values */
const std::string TYPE_ALLELES = "TYPE_ALLELES";			/* Mutation type for the alleles */
const std::string INIT_CLONAL = "INIT_CLONAL";	 			/* Initial clonal status */
const std::string ENVIRO_SDINIT = "ENVIRO_SDINIT";
const std::string ENVIRO_SDDYNAM = "ENVIRO_SDDYNAM";
const std::string ENVIRO_SDFINAL = "ENVIRO_SDFINAL";
const std::string ENVIRO_PLASTICITY = "ENVIRO_PLASTICITY";
const std::string FITNESS_TYPE = "FITNESS_TYPE";			/* Selection types */
const std::string FITNESS_STRENGTH = "FITNESS_STRENGTH";
const std::string FITNESS_OPTIMUM = "FITNESS_OPTIMUM";
const std::string FITNESS_CORRELATION = "FITNESS_CORRELATION";
const std::string FITNESS_STAB = "FITNESS_STAB";            /* Type of stability selection */
const std::string FITNESS_STABSTR = "FITNESS_STABSTR";      /* Strenght of selection on stability */
const std::string OUT_GENO = "OUT_GENO";					/* Output for the genotype */
const std::string OUT_UNSTAB = "OUT_UNSTAB";				/* Output for the phenotypic unstability */
const std::string OUT_CANAL_TESTS = "OUT_CANAL_TESTS";
const std::string OUT_CANAL_MUTSD = "OUT_CANAL_MUTSD";
const std::string OUT_CANAL_SDINIT = "OUT_CANAL_SDINIT";
const std::string OUT_CANAL_SDDYNAM = "OUT_CANAL_SDDYNAM";
const std::string OUT_HERIT_TESTS = "OUT_HERIT_TESTS";
const std::string OUT_DIREPI_TESTS= "OUT_DIREPI_TESTS";
const std::string PHENO_SCALING = "PHENO_SCALING";

// MULTILINEAR ARCHITECTURE
const std::string GENET_EPSILON2e = "GENET_EPSILON2e";
const std::string GENET_EPSILON2p = "GENET_EPSILON2p";
//const std::string GENET_EPSILON3 = "GENET_EPSILON3";

// BOOLEAN ARCHITECTURE
const std::string MATRIX_DENS = "MATRIX_DENS";
const std::string LOG_OPERATOR_DENS = "LOG_OPERATOR_DENS";
const std::string PHEN_NBLOC = "PHEN_NBLOC";
const std::string SCALE = "SCALE";

// REGULATORY ARCHITECTURE
const std::string INIT_CONNECT = "INIT_CONNECT";
const std::string INIT_CONDIAG = "INIT_CONDIAG";
const std::string TYPE_SO = "TYPE_SO";	
const std::string INIT_BASAL = "INIT_BASAL";				/* Initial vector (So) types */
const std::string INIT_RECURRENCE = "INIT_RECURRENCE";
const std::string DEV_TIMESTEPS = "DEV_TIMESTEPS";
const std::string DEV_CALCSTEPS = "DEV_CALCSTEPS";

// INPUT/OUTPUT
const std::string FILE_NEXTPAR = "FILE_NEXTPAR";


//Initial scale for boolean architecture
const std::string SC_int = "integer";
const std::string SC_vector = "vector";
const std::string SC_dec = "decimal";
const std::string SC_combi = "combined";
const std::vector<std::string> SC_options = boost::assign::list_of (SC_int)(SC_vector)(SC_dec)(SC_combi);

// Mutation type 
const std::string MT_individual = "individual";
const std::string MT_locus = "locus";
const std::vector<std::string> MT_options = boost::assign::list_of (MT_individual)(MT_locus);

// Mutation memory
const std::string MM_cumul = "cumul";
const std::string MM_stationary = "stationary";
const std::vector<std::string> MM_options = boost::assign::list_of (MM_cumul)(MM_stationary);

// Initial clonal status
const std::string CL_clonal = "clonal";
const std::string CL_variable = "variable";
const std::vector<std::string> CL_options = boost::assign::list_of (CL_clonal)(CL_variable);

// Selection types
const std::string FT_nosel = "no_selection";
const std::string FT_linear = "linear";
const std::string FT_expo = "exponential";
const std::string FT_gauss = "gaussian";
const std::string FT_multigauss = "multivar_gaussian";
const std::string FT_quad = "quadratic";
const std::string FT_truncup = "trunc_up";
const std::string FT_truncdown = "trunc_down";
const std::string FT_concave = "concave";
const std::string FT_convex = "convex";
const std::vector<std::string> FT_options = boost::assign::list_of (FT_nosel)(FT_linear)(FT_expo)(FT_gauss)(FT_multigauss)(FT_quad)(FT_truncup)(FT_truncdown)(FT_concave)(FT_convex);

// Strenght of selection on stability
const std::string FS_nostab = "no_stabsel";
const std::string FS_expo = "exponential_stab";
const std::vector<std::string> FS_options = boost::assign::list_of (FS_nostab)(FS_expo);

// Architecture types
const std::string AR_add = "additive";
const std::string AR_mult = "multilinear";
const std::string AR_wagner = "wagner";
const std::string AR_siegal = "siegal";
const std::string AR_m2 = "m2";
const std::string AR_Boolean = "boolean";
const std::vector<std::string> AR_options = boost::assign::list_of (AR_add)(AR_mult)(AR_wagner)(AR_siegal)(AR_m2)(AR_Boolean);

// Output for the genotype 
const std::string OG_yes = "yes";
const std::string OG_no = "no";
const std::vector<std::string> OG_options = boost::assign::list_of (OG_yes)(OG_no);

// Output for the phenotypic unstability
const std::string OU_yes = "yes";
const std::string OU_log = "log";
const std::string OU_no = "no";
const std::vector<std::string> OU_options = boost::assign::list_of (OU_yes)(OU_log)(OU_no);

// Initial vector type
const std::string SO_min = "minimum";
const std::string SO_max = "maximum";
const std::string SO_med = "median";
const std::string SO_randbin = "random_binary";
const std::string SO_rand = "random";
const std::string SO_basal = "basal";
const std::vector<std::string> SO_options = boost::assign::list_of (SO_min)(SO_max)(SO_med)(SO_randbin)(SO_rand)(SO_basal);

// Mutation type for the alleles
const std::string TA_norm = "normal";   // all mutations possible
const std::string TA_zero = "zero";     // if exactly 0.0: no mutations, otherwise mutations possible
const std::string TA_immut= "immut";    // the allele can never mutate and stays at its initial state
const std::string TA_sign = "sign";     // mutations are possible, but the sign is unchanged
const std::vector<std::string> TA_options = boost::assign::list_of (TA_norm)(TA_zero)(TA_immut)(TA_sign);

// Scaling transformations for the phenotype
const std::string ST_none = "none";
const std::string ST_log  = "log";
const std::string ST_logit= "logit";
const std::string ST_m1101= "m1101";
const std::string ST_invlogit="invlogit";
const std::vector<std::string> ST_options = boost::assign::list_of (ST_none)(ST_log)(ST_logit)(ST_m1101)(ST_invlogit);


#endif // PARCONST_H_INCLUDED
