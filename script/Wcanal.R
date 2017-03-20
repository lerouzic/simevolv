#########################################################
# Wcanal.R
#
# Computes the canalization properties of a W matrix
# 
# Copyright Arnaud Le Rouzic / CNRS 2015
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


script.dir <- normalizePath(dirname(parent.frame(2)$ofile))
if (is.null(script.dir))
	script.dir <- normalizePath(dirname(sys.frame(1)$ofile))
	
source(paste(script.dir, "netw.R", sep="/"))
source(paste(script.dir, "misc.R", sep="/"))


dev.stability <- function(W, ...) {
	# Internal (genetic) stability at the end of the development
	canIndex(model.M2(W, ...)$var)
}

homeostasis <- function(W, a=0.5, S0=rep(a, nrow(W)), ..., microSdE=0.1, nb.affected.genes=ncol(W), delay=1, replicates=100) {
	# Influence of Gaussian microenvironmental perturbations at the end of the development
	# on the mean expression
    # In practice, the score is the square of the difference between 
    ref <- model.M2(W, a=a, S0=S0, ...)$mean
	rr <- mcsapply(1:replicates, function(i) {
        newS0 <- ref
        affected.genes <- sample.int(length(newS0), nb.affected.genes)
        newS0[affected.genes] <- newS0[affected.genes] + rnorm(nb.affected.genes, 0, sd=microSdE)
        newS0[newS0 < 0] <- 0
        newS0[newS0 > 1] <- 1
        model.M2(W, a=a, S0=newS0, steps=delay, measure=1, ...)$mean })
    # Computes the mean square of the differences before/after disturbance
	canIndex(apply((rr-ref)^2, 1, mean))
}

dev.canalization <- function(W, a=0.5, S0=rep(a, nrow(W)), ..., nb.affected.genes=ncol(W), replicates=100) {
	# Influence of the initial expression on the final mean expression

    ref <- model.M2(W, a=a, S0=S0, ...)$mean
	rr <- mcsapply(1:replicates, function(i) {
        newS0 <- ref
        affected.genes <- sample.int(length(newS0), nb.affected.genes)
        newS0[affected.genes] <- runif(nb.affected.genes)
		model.M2(W=W, a=a, S0=newS0, ...)$mean
		})
    # Computes the mean square of the differences before/after disturbance
	canIndex(apply((rr-ref)^2, 1, mean))
}

canscale <- function(x) ifelse(is.na(x) | x < exp(-20), 20, -log(x))

genet.canalization <- function(W, ..., which.i=1:nrow(W), which.j=1:ncol(W), mutsd=0.1, exclude.0=TRUE, replicates=100, what="mean", reference="mutant") {
	# This is the traditional canalization measurement (robustness to mutations)
	# reference = 'mutant' (variance around the average mutant expression)
	# reference = 'wt'     (variance around the original genotype)
	#
	# which.i and which.j specify the lines/columns that can be mutated
	# exclude.0 : disallow mutations of the 0s
	
	stopifnot(nrow(W) == ncol(W))
	n <- nrow(W)
	stopifnot(all(which.i > 0 & which.i <= n))
	stopifnot(all(which.j > 0 & which.j <= n))	
	
	mrow <- matrix(rep(1:n, n), nrow=n)
	mcol <- t(mrow)
	
	avoided <- !(mrow%in%which.i) | !(mcol%in%which.j) | (W == 0 & exclude.0)
	stopifnot (sum(avoided) < length(W))

	Wbox <- if (sum(!avoided) == 1) { rep(which(!avoided), replicates)}
			else {sample(which(!avoided), replicates, replace=TRUE)}
	Wdev <- rnorm(replicates, 0, sd=mutsd)
	rr <- mcsapply(1:replicates, function(i) {
		myW <- W
		W[Wbox[i]] <- W[Wbox[i]]+Wdev[i]
		model.M2(W=W, ...)$mean
	})

	if (reference == "mutant") {
		return(apply(rr, 1, var))
	} else if (reference == "wt") {
		ref <- model.M2(W=W, ...)[[what]]
		return(apply(rr-ref, 1, function(x) sum(x^2)/length(x)))
	}	
}

