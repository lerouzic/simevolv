#!/usr/bin/env Rscript
#
#################################################################################################
# R script used to generate parameter files for evolutive simulations with "simevolv".		#
# These parameter files only contain an optimum and the link to the next parameter file		#
# This is used to produce the series of optima to be reached in flcutuating selection		#
# 												#
# FURTHER versions of this script should allow the user to change :				#
# 	distribution and range  ;  signal-opt correlation  ;  modification to the signal	#
#		modifications to the signal including : constant/random signal or noisy signal	#
# 												#
# 												#
# Copyright Andreas Odorico / Paris-Sud University 2018						#
# <andreas.odorico@u-psud.fr>									#
# 												#
# Released under the WTFPL version 2.0								#
# * No warranty *										# 
# 												#
#################################################################################################

## LIST OF ARGUMENTS: (See "help" for additional informations)

# REPLICATES (-rep) ; GENERATIONS (-gen) ; envSignal DISTRIBUTION (and range) (-distrib) ; CORRELATION (signal-opt) (-corr)
# DIGITS (-digits) ; signalModif_Type (-smType) ; signalModif_Value (-smVal) ; signalNOISE (-smNoise)

##################################################################################################

## Global Warning :

### Inputs are not (yet) checked. Be careful.

### Script is hardwire for only ONE signal gene which is the FIRST

##################################################################################################
## Requirements.
require("parallel")

## set Defaut Parameters

NB_CORES <- 1
#PAATH		#no default value
#nREP		#no default value
#nGENER		#no default value
namePTRN <- 'paste0("par_rep",rrr,"_gen",ggg)'

DISTRIB <- "runif(1,0,1)"
DIGITS <- 3
CORRELATION <- "x ; x ; 1-x ; ((1/2)+(x/2)) ; ((1/2)-(x/2)) ; 0.5 ; 0.5 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0"

signalModif_Type <- FALSE
	# FALSE => signal has no difference with optimum.
	# CST	=> overwrite signal value with "signalModif_Value"
	# RNDM	=> will pick signal in the distribution given in "signalModif_Value"
signalModif_Value <- runif(1,0,1)	# would be CST (same for ALL SIMULATIONS, ALL GENERATIONS) # I may want a diff btw sims...
signalModif_Value <- "runif(1,0,1)"	# would be RNDM (uniform between 0 et 1 for EACH GENERATION, same for all pops)

signalNOISE <- NULL
#signalNOISE <- "rnorm(1,0,0.2)" # expected
	# I could just set a double which would be the SD of the rnorm... as no other distribution really makes sense here
##################################################################################################
## Read Arguments :

if(!exists("ARGUMENTS")) { ARGUMENTS <- commandArgs(trailingOnly=TRUE) }
	if(length(which(ARGUMENTS=="-proc"))!=0) { NB_CORES <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-proc")+1]) }
	if(length(which(ARGUMENTS=="-path"))!=0) { PAATH <- ARGUMENTS[which(ARGUMENTS=="-path")+1] }
	if(length(which(ARGUMENTS=="-rep"))!=0) { nREP <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-rep")+1]) }
	if(length(which(ARGUMENTS=="-gen"))!=0) { nGENER <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-gen")+1]) }

	if(length(which(ARGUMENTS=="-name"))!=0) { namePTRN <- ARGUMENTS[which(ARGUMENTS=="-name")+1] }
	if(length(which(ARGUMENTS=="-distrib"))!=0) { DISTRIB <- ARGUMENTS[which(ARGUMENTS=="-distrib")+1] }
	if(length(which(ARGUMENTS=="-digits"))!=0) { DIGITS <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-digits")+1]) }
	if(length(which(ARGUMENTS=="-corr"))!=0) { CORRELATION <- ARGUMENTS[which(ARGUMENTS=="-corr")+1] }
	if(length(which(ARGUMENTS=="-smType"))!=0) { signalModif_TYPE <- ARGUMENTS[which(ARGUMENTS=="-smType")+1] }
	if(length(which(ARGUMENTS=="-smVal"))!=0) { signalModif_VALUE <- ARGUMENTS[which(ARGUMENTS=="-smVal")+1] }
	if(length(which(ARGUMENTS=="-smNoise"))!=0) { signalNOISE <- ARGUMENTS[which(ARGUMENTS=="-smNoise")+1] }


#(test)
#ARGUMENTS <- c("-rep", "12", "-proc", "75", "-path", "./C://mescouilles/", "-distrib", "rpoiss(0,666,1337)", "-gen", 2, "-gene", "-argument", "--test", "shama-graba", "-howbout-this?")

called_args <- ARGUMENTS[grep("^-", ARGUMENTS)]
existing_args <- c("-proc", "-path", "-rep", "-gen", "-name", "-distrib", "-digits", "-corr", "-smType", "-smVal", "-smNoise", "-h", "--help")

if(sum(!(called_args %in% existing_args))>0) {
	uncorrect_args <- c("{ ", paste(called_args[!(called_args %in% existing_args)], collapse=" ; "), " }")
	stop(c("Argument(s) : ", uncorrect_args, " is (are) not recognized. Use \"-h\" or \"--help\" to display a full list of recognized arguments along with usage informations."))
}


### Help message :
if( (length(ARGUMENTS[which(ARGUMENTS=="-h")])!=0) || (length(ARGUMENTS[which(ARGUMENTS=="-help")])!=0) ) {
	message(c("\n\tThis script generates linked parameter files with flucutating optima in them. All fluctuating optima are calculated as functions of an environmental signal. By default, this signal is also used to compute the expression value of a signal gene that serves as input for a plastic regulation network.\n",
"Parameter files will be produced for each \"-rep\" population and \"-gen\" generation, following the \"-name\" pattern.\n\n",
"Only \"-path\", \"-rep\" and \"-gen\" have no default values. Please PROVIDE ALL ARGUMENTS AS STRINGS INTERPRETABLE BY R \"eval(parse())\" and WITHOUT SPACES for a given argument.\n",

"This script only supports using the FIRST GENE AS A SIGNAL.\n\n",
"Also few assertions are made... so please be careful when using exotic arguments. \nSupported Arguments :\n",

"\t\"-path\" must be provided. Defines where paramFiles must be created.\n\n",
"\t\"-rep\" must be provided. It sets the number of populations for which optima must be generated\n\n",
"\t\"-gen\" must be provided. It sets the number of generations for which optima must be generated\n\n",
"\t\"-name\" sets the pattern of the name for the parameter files created. By default, \"par_repX_genY\" where \"X\" is replicate number and \"Y\" is generation number.\n\t\tIn the R code, \"X\" is \"rrr\" and \"Y\" is \"ggg\" so the default value is the string \"paste0(\"par_rep\",rrr,\"_gen\",ggg)\". Again, without spaces.\n\t\tPlease also note that the \"ggg\" iterator must be recognized (with whole word matching) to set the next_generation_parfile_name. Avoid any ingenious tricks to confuse \"grep(\\\\bggg\\\\b)\".\n\n",

"\t\"-distrib\" sets the distribution and range of environmental signal values used to compute optima and signals. Default value is \"runif(1,0,1)\".\n\n",
"\t\"-digits\" sets the precision for the optimum values and environmental signal. Default is 3.\n\n",
"\t\"-corr\" defines the correlation between optima and environmental signal. It must be provided as a series of functions of \"x\" separated by SEMICOLONS (\";\").\n\t\tDefault is \"x ; x ; 1-x ; ((1/2)+(x/2)) ; ((1/2)-(x/2)) ; 0.5 ; 0.5\" followed by 13 zeroes.\n\n",
"\t\"-smType\" is FALSE by default, \"CST\" or \"RNDM\" overwrite the signal with the value given by the argument \"-smVal\"\n\n",
"\t\"-smVal\" : if \"-smType\" is not False, requires a number for \"CST\" signal overwriting or a distribution for \"RNDM\" signal generation.\n\n",
"\t\"-smNoise\" is NULL by default. May be filled with anything that can be evaluated as a single value. The expression of the signal gene (1) will be modified by this value.\n\n",
"\t\"-proc\" sets the number of processing cores used. Uses R \"parallel\" package. Default is 1. Removing the dependency on the \"parallel\" package even for NB_CORES=1 requires a bit of recoding...\n"))
	quit("no")	
}

if(!exists("nREP") || length(nREP)==0) stop("Argument -rep not supplied and has no default value")
if(!exists("nGENER") || length(nGENER)==0) stop("Argument -gen not supplied and has no default value")
if(!exists("PAATH") || length(PAATH)==0) stop("Argument -path not supplied and has no default value")

##################################################################################################
## Script :

invisible(lapply(1:nREP, function(rrr) {
	invisible(lapply(2:nGENER, function(ggg) {	
# pick an environment-signal in the provided distribution
	envSignal <- eval(parse(text=DISTRIB))

# Define all optima as the provided "correlation function" of environment-signal
	OPTIMA <- unlist(strsplit(x=CORRELATION, split=";"))
	OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(envSignal, function(x) eval(parse(text=OPTIMA[oo])) ) )

# Overwrite signal-gene optimum (by default, set to environment-signal)
	if(signalModif_Type) {
		if(signalModif_Type=="CST") { OPTIMA[1] <- as.numeric(eval(parse(text=signalModif_Value)))
		} else if(signalModif_Type=="RNDM") { OPTIMA[1] <- as.numeric(eval(parse(text=signalModif_Value))) }
	}
# Add noise to signal-gene value:
	if(!is.null(signalNOISE)) { OPTIMA[1] <- OPTIMA[1] + as.numeric(eval(parse(text=signalNOISE))) }

# Writing result:
	OPTIMA <- round(OPTIMA, digits=DIGITS)
		#some assertions:
	OPTIMA <- as.numeric(OPTIMA)
	if( anyNA(OPTIMA) || length(which(OPTIMA>1))>1 || length(which(OPTIMA<0))>1 ) { stop("At least one optimum is either \"NA\" or out of the 0-1 range. Check \"-distrib\" or \"-corr\"") }

	ANS <- c("FITNESS_OPTIMUM", paste(OPTIMA, collapse=" "))
	#Last generation : no FILE_NEXTPAR
	if(ggg < nGENER) {
		NEXT_PARFILE <- eval(parse(text= gsub("\\bggg\\b","(ggg+1)",namePTRN) ))
		ANS <- c(ANS, "FILE_NEXTPAR", NEXT_PARFILE)
	}
	write(ANS, file=paste(PAATH, eval(parse(text=namePTRN)), sep="/"), sep="\t\t", ncolumns=2)
	})) #Generation Loop
})) #Replicates Loop


