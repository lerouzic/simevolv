#########################################################
# readfile.R
#
# Reads simulation files and extract information
# (such as W matrices)
# 
# Copyright 
#    Arnaud Le Rouzic 
#    Estelle Runneburger
#    Andreas Odorico  
#    CNRS-Univ. Paris-Sud 2015-2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


read.W <- function(simfile, gen=NA, last.only=FALSE, first=TRUE, pattern="MeanAll") {
	stopifnot(
		is.character(simfile),
		file.exists(simfile)
	)
	
	daata <- read.table(simfile, header=TRUE)
	if (is.na(gen)) {
		gen <- daata$Gen
		if (!first) gen <- gen[-1]
		if (last.only) gen <- gen[length(gen)]
	}
	ans <- list()
	for (gg in gen) {
		mat <- as.numeric(daata[which(daata$Gen == gg), grep(pattern, colnames(daata))])
		if (length(mat) == 0) stop("No field called \"", pattern, "\" in file \"", simfile, "\".")
		mat <- matrix(mat, ncol=sqrt(length(mat)), byrow=TRUE) #byrow=FALSE to have t(W)
		ans[[as.character(gg)]] <- mat
	}
#~	return(ans)
	return(lapply(ans, function(a) unlist(a, recursive=FALSE, use.names=FALSE)))
}


read.W_csv <- function(simfile, gen=NA, last.only=FALSE, first=TRUE, pattern="MeanAll") {
	stopifnot(
		is.character(simfile),
		file.exists(simfile)
	)
	daata <- read.csv(simfile, header=TRUE)
	if (is.na(gen)) {
		gen <- daata$Gen
		if (!first) gen <- gen[-1]
		if (last.only) gen <- gen[length(gen)]
	}
	ans <- list()
	for (gg in gen) {
		mat <- as.numeric(daata[which(daata$Gen == gg), grep(pattern, colnames(daata))])
		if (length(mat) == 0) stop("No field called \"", pattern, "\" in file \"", simfile, "\".")
		mat <- matrix(mat, ncol=sqrt(length(mat)), byrow=TRUE) #byrow=FALSE to have t(W)
		ans[[as.character(gg)]] <- mat
	}
	return(ans)
}

