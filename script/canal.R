# #########################################################
# canal.R
#
# Various ways to analyse the evolution of canalization scores
# 
# Copyright Arnaud Le Rouzic / CNRS 2016
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
#
#
# Quick and Dirty documentation:
# 
# example of usage in R:
# 
# DeltaCan(c("file1.txt", "file2.txt", "file3.txt"), from=5000)
#
# options:
# * method = "lm" or "raw". "lm" performs a linear regression, while "raw" computes a difference.
# * from = generation from which the DeltaC is calculated. 
# * outcan = "can" or "mut" according to the option in the simulation parameter
# * avg.files = whether the deltaC should be averaged over the files
# * avg.genes = whether the deltaC should be averaged over the genes
#
# In case "from" is missing, an educated guess is made. The calculation starts from the timepoint
# following the minimal canalization. This has a number of drawbacks, mainly
# - if the time series is decreasing, returns NA
# - if several files/genes are considered, the starting point may be different and thus
#   the results will not be completely comparable
# - using the timepoint after the minimum decreases the bias due to starting from the lowest point, 
#   but this assumes no correlation between consecutive canalization measurements. 
#############################################################


DeltaCan <- function(files, method="lm", from=NA, outcan="can", avg.files=TRUE, avg.genes=FALSE) 
{
	dd <- lapply(files, function(ff) {
		dat <- read.table(ff, header=TRUE)
		dat[,c(1, grep("MCan", colnames(dat)))] })
	
	if (outcan == "mut")
		dd <- lapply(dd, function(d) cbind(d[,1], -log(d[,-1])))

	cs <- lapply(dd, canscore, method=method, from=from)
	cs <- do.call(rbind, cs)
	rownames(cs) <- files
	colnames(cs) <- paste0("MCan", 1:ncol(cs))
	if (avg.genes) cs <- apply(cs, 1, mean)
	if (avg.files) if (is.matrix(cs)) cs <- apply(cs, 2, mean) else cs <- mean(cs)
	return(cs)
}

canscore <- function(x, method, from) {
	if (method=="lm") {
		return(sapply(2:ncol(x), function(i) canscore.lm(x[,i], from, gen=x[,1])))
	} else if (method=="raw") {
		return(sapply(2:ncol(x), function(i) canscore.raw(x[,i], from, gen=x[,1])))
	} else {
		stop("Method", method, "unknown.")
	}
}

canscore.raw <- function(x, from, gen) {
	gen <- as.numeric(gen)
	stopifnot(
		length(x) == length(gen),
		is.na(from) || (from >= min(gen) && from < max(gen) && from %in% gen))
	if (is.na(from)) {
		# Guessing: from the timestep following the minimum.
		wm <- which.min(x)
		if (wm > length(x)-2) {
			warning("Decreasing time series, impossible to guess the beginning.")
			return(NA)
		}
		from <- gen[wm+1]
	}
	fm <- which(gen==from)
	return((x[length(x)] - x[fm])/(gen[length(x)] - from))
}

canscore.lm <- function(x, from, gen) {
	gen <- as.numeric(gen)
	stopifnot(
		length(x) == length(gen),
		is.na(from) || (from >= min(gen) && from < max(gen) && from %in% gen))
	if (is.na(from)) {
		# Guessing: from the timestep following the minimum.
		wm <- which.min(x)
		if (wm > length(x)-2) {
			warning("Increasing time series, impossible to guess the beginning.")
			return(NA)
		}
		from <- gen[wm+1]
	}
	ll <- lm(x[from <= gen] ~ gen[from <= gen])
	return(coef(ll)[2])
}
