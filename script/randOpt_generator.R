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
# DIGITS (-digits) ; signalModif_TYPE (-smType) ; signalModif_Value (-smVal) ; signalNOISE (-smNoise)

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
namePTRN <- 'paste0("par_rep",rrr,"_gen",ggg)'

DISTRIB <- "runif(1,0,1)"
DIGITS <- 3
CORRELATION <- "x ; x ; 1-x ; ((1/2)+(x/2)) ; ((1/2)-(x/2)) ; 0.5 ; 0.5 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0"

signalModif_TYPE <- "None"
	# CST	=> overwrite signal value with "signalModif_VALUE"
	# RNDM	=> will pick signal in the distribution given in "signalModif_VALUE"
	# any other value will be ignored.
signalModif_VALUE <- "runif(1,0,1)"
	# if signalModif_TYPE is CST and signalModif_VALUE is a constant, a same value will be used for all replicates.
	# an expression may be used to compute a given value once per replicate (with "CST") or once per generation ("RNDM")

signalNOISE <- NULL
#signalNOISE <- "rnorm(1,0,0.2)" # example.
	# I could just set a double which would be the SD of the rnorm... as no other distribution really makes sense here
##################################################################################################
## Read Arguments :

if(!exists("ARGUMENTS")) { ARGUMENTS <- commandArgs(trailingOnly=TRUE) }	
	if(length(which(ARGUMENTS=="-proc"))!=0) { NB_CORES <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-proc")+1]) }

	if(length(which(ARGUMENTS=="-param"))!=0) { G1PARAM <- ARGUMENTS[which(ARGUMENTS=="-param")+1] }
	if(length(which(ARGUMENTS=="-path"))!=0) { PAATH <- ARGUMENTS[which(ARGUMENTS=="-path")+1] }
	if(length(which(ARGUMENTS=="-rep"))!=0) { nREP <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-rep")+1]) }

	if(length(which(ARGUMENTS=="-name"))!=0) { namePTRN <- ARGUMENTS[which(ARGUMENTS=="-name")+1] }
	if(length(which(ARGUMENTS=="-distrib"))!=0) { DISTRIB <- ARGUMENTS[which(ARGUMENTS=="-distrib")+1] }
	if(length(which(ARGUMENTS=="-digits"))!=0) { DIGITS <- as.numeric(ARGUMENTS[which(ARGUMENTS=="-digits")+1]) }
	if(length(which(ARGUMENTS=="-corr"))!=0) { CORRELATION <- ARGUMENTS[which(ARGUMENTS=="-corr")+1] }
	if(length(which(ARGUMENTS=="-smType"))!=0) { signalModif_TYPE <- ARGUMENTS[which(ARGUMENTS=="-smType")+1] }
	if(length(which(ARGUMENTS=="-smVal"))!=0) { signalModif_VALUE <- ARGUMENTS[which(ARGUMENTS=="-smVal")+1] }
	if(length(which(ARGUMENTS=="-smNoise"))!=0) { signalNOISE <- ARGUMENTS[which(ARGUMENTS=="-smNoise")+1] }


called_args <- ARGUMENTS[grep("^-", ARGUMENTS)]
existing_args <- c("-param", "-path", "-rep", "-proc", "-name", "-distrib", "-digits", "-corr", "-smType", "-smVal", "-smNoise", "-h", "--help")

if(sum(!(called_args %in% existing_args))>0) {
	uncorrect_args <- c("{ ", paste(called_args[!(called_args %in% existing_args)], collapse=" ; "), " }")
	stop(c("Argument(s) : ", uncorrect_args, " is (are) not recognized. Use \"-h\" or \"--help\" to display a full list of recognized arguments along with usage informations."))
}

### Help message :
if( (length(ARGUMENTS[which(ARGUMENTS=="-h")])!=0) || (length(ARGUMENTS[which(ARGUMENTS=="--help")])!=0) ) {
	message(c("\n\tThis script generates linked parameter files with flucutating optima in them. All fluctuating optima are calculated as functions of an environmental signal. By default, this signal is also used to compute the expression value of a signal gene that serves as input for a plastic regulation network.\n",
"Parameter files will be produced for each \"-rep\" population, following the \"-name\" pattern. The length and period of fluctuations are deduced from the parameter file (respectively SIMUL_MAXGEN and SIMUL_GENER). With this generator, the period of fluctuations IS FIXED. For fluctuations of heterogeneous length, please devise another script or repeat this one with saving and loading of populations between every set of parameter files.\n\n",
"Only \"-path\", \"-rep\" and \"-param\" have no default values. Please PROVIDE ALL ARGUMENTS AS STRINGS INTERPRETABLE BY R \"eval(parse())\" and WITHOUT SPACES for a given argument.\n",

"This script only supports using the FIRST GENE AS A SIGNAL.\n\n",
"Also few assertions are made... so please be careful when using exotic arguments. \nSupported Arguments :\n",

"\t\"-param\" must be provided. It mainly used to set the \"generation 1\" parameter file. Which is a copy of the provided file with modified optimum values to match the \"distribution\" and \"correlation\" given, as are the following generations. Also a \"FILE_NEXTPAR\" argument is added in line one, make sure no other value of FILE_NEXTPAR is already provided in the parameter file.\n Also, the \"SIMUL_GENER\" parameter in the parameter file will be used to set the period of fluctuations.\n\n",

"\t\"-path\" must be provided. Preferably Absolute Path. Defines where paramFiles must be created.\n\n",
"\t\"-rep\" must be provided. It sets the number of populations for which optima must be generated\n\n",
"\t\"-name\" sets the pattern of the name for the parameter files created. By default, \"par_repX_genY\" where \"X\" is replicate number and \"Y\" is generation number.\n\t\tIn the R code, \"X\" is \"rrr\" and \"Y\" is \"ggg\" so the default value is the string \"paste0(\"par_rep\",rrr,\"_gen\",ggg)\". Again, without spaces.\n\t\tPlease also note that the \"ggg\" generation index must be recognized (with whole word matching) to set the next_generation_parfile_name. Avoid any ingenious tricks to confuse \"grep(\\\\bggg\\\\b)\" (the latter being used to .\n\n",

"\t\"-distrib\" sets the distribution and range of environmental signal values used to compute optima and signals. Default value is \"runif(1,0,1)\".\n\n",
"\t\"-digits\" sets the precision for the optimum values and environmental signal. Default is 3.\n\n",
"\t\"-corr\" defines the correlation between optima and environmental signal. It must be provided as a series of functions of \"x\" separated by SEMICOLONS (\";\").\n\t\tDefault is \"x ; x ; 1-x ; ((1/2)+(x/2)) ; ((1/2)-(x/2)) ; 0.5 ; 0.5\" followed by 13 zeroes.\n\n",
"\t\"-smType\" \"CST\" or \"RNDM\" overwrite the signal value with the value given by the argument \"-smVal\". Any other value will be ignored (default is \"None\").\n\n",
"\t\"-smVal\" : if \"-smType\" is not \"NONE\", requires a distribution for \"RNDM\" signal generation and to generate a different \"CST\" signal for every replicate. If a number is provided and \"-smType\" is \"CST\", the same signal will be used for all replicates.\n\n",
"\t\"-smNoise\" is NULL by default. May be filled with anything that can be evaluated as a single value. The expression of the signal gene (1) will be modified by this value.\n\n",
"\t\"-proc\" sets the number of processing cores used. Uses R \"parallel\" package. Default is 1. Removing the dependency on the \"parallel\" package even for NB_CORES=1 requires a bit of recoding...\n"))
	quit("no")	
}

if(!exists("nREP") || length(nREP)==0) stop("Argument -rep not supplied and has no default value")
if(!exists("PAATH") || length(PAATH)==0) stop("Argument -path not supplied and has no default value")
if(!exists("G1PARAM") || length(G1PARAM)==0) stop("Argument -param not supplied and has no default value")

if(signalModif_TYPE=="RNDM" && is.numeric(signalModif_VALUE)) stop("\"-smType\" is \"RNDM\", please use a string containing a distribution interpretable by R for \"-smVal\". Or set \"-smType\" to \"CST\".")

##################################################################################################
## Script :
### Assess at which generations fluctuations are to be made.
cmd <- paste0("echo ","$(awk \'/SIMUL_MAXGEN/{print $NF}' ", G1PARAM,")")
	GENMAX <- as.numeric(system(cmd, intern=TRUE))
cmd <- paste0("echo ","$(awk \'/SIMUL_GENER/{print $NF}' ", G1PARAM,")")
	GENSTEP<- as.numeric(system(cmd, intern=TRUE))
GENSEQ <- seq(0, GENMAX, GENSTEP)
if(GENSEQ[length(GENSEQ)]!=GENMAX) { GENSEQ <- append(GENSEQ, GENMAX) }
GENSEQ <- GENSEQ[-1]				# Optimum parameterFiles created are incomplete. The first one need to be a copy of a complete parameter file.
if(GENSEQ[1]==1) { GENSEQ <- GENSEQ[-1] }
	## Sequence of generations goes : "1, 1+N", 1+(2N), 1+(3N)..." => as the "1" is fixed, it is removed from GENSEQ on which "ggg" will iterate AFTER generation "1"
		# the "step" of the sequence of generations is fixed with this generator. i.e. "1, 1+N, 1+N+M, 1+N+M+Y..." cannot be generated with this script in a single run, as is.

### Generate linked optimum files
invisible(lapply(1:nREP, function(rrr) {
	# copy the provided parameterfile ; run once the process of optimum generation, write it ; add FILE_NEXTPAR on the first line.
		envSignal <- eval(parse(text=DISTRIB))
		OPTIMA <- unlist(strsplit(x=CORRELATION, split=";"))
		OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(envSignal, function(x) eval(parse(text=OPTIMA[oo])) ) )

		if(signalModif_TYPE=="CST") { TEMP_CST_signVal <- as.numeric(eval(parse(text=signalModif_VALUE))) ; OPTIMA[1] <- TEMP_CST_signVal
		} else if(signalModif_TYPE=="RNDM") { OPTIMA[1] <- as.numeric(eval(parse(text=signalModif_VALUE)))	# "RDNM" will evaluate the distribution at every generation
		} else {} #any other value does not overwrite the signal... This line is useless, but explicit.

		if(!is.null(signalNOISE)) { OPTIMA[1] <- OPTIMA[1] + as.numeric(eval(parse(text=signalNOISE))) }
	# Writing
		OPTIMA <- round(OPTIMA, digits=DIGITS)
		OPTIMA <- as.numeric(OPTIMA)
			OPTIMA[which(OPTIMA>1)] <- 1
			OPTIMA[which(OPTIMA<0)] <- 0
		if( anyNA(OPTIMA) | length(which(OPTIMA>1))>0 | length(which(OPTIMA<0))>0 ) { stop("At least one optimum is either \"NA\" or out of the 0-1 range. Check \"-distrib\" or \"-corr\"") }
	OPTIMA <- paste(OPTIMA, collapse=" ")
	# copy the param file
	cmd <- paste0("cp ",G1PARAM," ",paste(PAATH, eval(parse(text= gsub("\\bggg\\b","1",namePTRN) )), sep="/"))
	system(cmd)
	# replace FITNESS_OPTIMUM in first line
	cmd <- paste0("sed -i -e \"s/FITNESS_OPTIMUM[[:blank:]].*/FITNESS_OPTIMUM\t",OPTIMA,"/\" ", paste(PAATH, eval(parse(text= gsub("\\bggg\\b","1",namePTRN) )), sep="/"))
	system(cmd)
	# add FILE_NEXTPAR in first line (didn't even try to do it in one file access)
	cmd <- paste0("sed -i \"1 i\\FILE_NEXTPAR\t",paste(PAATH, eval(parse(text= gsub("\\bggg\\b","GENSEQ[1]",namePTRN) )), sep="/"),"\" ",
			paste(PAATH, eval(parse(text= gsub("\\bggg\\b","1",namePTRN) )), sep="/") )
	system(cmd)

######## LOOP OVER THE GENERATIONS
	mclapply(GENSEQ, function(ggg) {
# pick an environment-signal in the provided distribution
	envSignal <- eval(parse(text=DISTRIB))

# Define all optima as the provided "correlation function" of environment-signal
	OPTIMA <- unlist(strsplit(x=CORRELATION, split=";"))
	OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(envSignal, function(x) eval(parse(text=OPTIMA[oo])) ) )

# Overwrite signal-gene optimum (by default, set to environment-signal)
		if(signalModif_TYPE=="CST") { OPTIMA[1] <- TEMP_CST_signVal
	} else if(signalModif_TYPE=="RNDM") { OPTIMA[1] <- as.numeric(eval(parse(text=signalModif_VALUE)))
	} else {} # Again, useless line, kept for explicity.

# Add noise to signal-gene value:
	if(!is.null(signalNOISE)) { OPTIMA[1] <- OPTIMA[1] + as.numeric(eval(parse(text=signalNOISE))) }

# Writing result:
	OPTIMA <- round(OPTIMA, digits=DIGITS)
		#some assertions:
	OPTIMA <- as.numeric(OPTIMA)
			OPTIMA[which(OPTIMA>1)] <- 1
			OPTIMA[which(OPTIMA<0)] <- 0
	if( anyNA(OPTIMA) | length(which(OPTIMA>1))>0 | length(which(OPTIMA<0))>0 ) { stop("At least one optimum is either \"NA\" or out of the 0-1 range. Check \"-distrib\" or \"-corr\"") }

	ANS <- c("FITNESS_OPTIMUM", paste(OPTIMA, collapse=" "))
	#Last generation : no FILE_NEXTPAR
	if(ggg < GENMAX) {
		NEXT_PARFILE <- eval(parse(text= gsub("\\bggg\\b","(GENSEQ[which(GENSEQ==ggg)+1])",namePTRN) ))
		ANS <- c(ANS, "FILE_NEXTPAR", paste(PAATH, NEXT_PARFILE, sep="/") )
	}
	write(ANS, file=paste(PAATH, eval(parse(text=namePTRN)), sep="/"), sep="\t\t", ncolumns=2)

	}, mc.cores=NB_CORES) #Generation Loop
})) #Replicates Loop

