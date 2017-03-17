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

dev.stability <- function(W, ...) {
	# Internal (genetic) stability at the end of the development
	model.M2(W, ...)$var
}

dev.robustness <- function(W, ..., microSdE=0.1, replicates=100, what="mean") {
	# Influence of microenvironmental perturbations during development
	# on the mean expression
	rr <- sapply(1:replicates, function(i) model.M2(W, ..., microSdE=microSdE)[[what]])
	apply(rr, 1, var)
}

dev.consistency <- function(W, a=0.5, S0=rep(a, nrow(W)), ..., macroSdE=0.1, replicates=100, what="mean") {
	# Influence of the initial expression on the final mean expression

	rr <- sapply(1:replicates, function(i) {
		S <- rnorm(length(S0), mean=S0, sd=macroSdE)
		S[S < 0] <- 0
		S[S > 1] <- 1		
		model.M2(W=W, a=a, S0=S, ...)[[what]]
		})
	apply(rr, 1, var)
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
	rr <- sapply(1:replicates, function(i) {
		myW <- W
		myW[Wbox[i]] <- W[Wbox[i]]+Wdev[i]
		model.M2(W=myW, ...)[[what]]
	})
	if (reference == "mutant") {
		return(apply(rr, 1, var))
	} else if (reference == "wt") {
		ref <- model.M2(W=W, ...)[[what]]
		return(apply(rr-ref, 1, function(x) sum(x^2)/length(x)))
	}	
}

genet.canalization.ij <- function(W, ..., mutsd=0.1, replicates=100, what="mean", reference="mutant", plot=FALSE, on.gene=1:nrow(W), zlim=c(-log(0.25), 20), text=TRUE) {
	# Computes the canalization score for every single regulatory element (i, j)
	dat <- expand.grid(1:nrow(W), 1:nrow(W))
	colnames(dat) <- c("i","j")
	dat <- cbind(dat, t(mapply(dat$i, dat$j, FUN=function(i, j) 
		canscale(genet.canalization(W=W, ..., which.i=i, which.j=j, 
			exclude.0=FALSE, mutsd=mutsd, replicates=replicates, 
			reference=reference, what=what)))
		)
	)
	if (plot) {
		layout(t(1:2), width=c(0.8, 0.2))
		mcan <- matrix(apply(dat[,-c(1,2)], 1, function(x) mean(x[on.gene])), ncol=ncol(W), byrow=FALSE)
		n <- nrow(mcan)
		image(t(mcan)[,n:1], x=1:n, y=1:n, zlim=zlim, xaxt="n", yaxt="n", col=heat.colors(100))
		axis(1, at=1:n, as.character(1:n))
		axis(2, at=n:1, as.character(1:n))
		if(text) text(x=dat$j, y=n+1-dat$i, labels=as.character(round(mcan, digits=1)))
		
		ref <- model.M2(W=W, ...)[[what]]
		
		image(t(as.matrix(rev(ref))), y=1:n, xaxt="n", yaxt="n", col=colorRampPalette(c("white","blue"))(100))
		if (text) text(x=0, y=n:1, labels=round(ref, 2))
		layout(1)
	}
	return(dat)
}

test.canalization <- function(W, a=0.5, S0=rep(a, nrow(W)), ..., microSdE=0.1, macroSdE=0.1, mutsd=0.1, replicates=100) {
	ans <- cbind(
		Stability=dev.stability(W, ...), 
		Robustness=dev.robustness(W, ..., microSdE=microSdE, replicates=replicates),
		Consistency=dev.consistency(W, a=a, S0=S0, ..., macroSdE=macroSdE, replicates=replicates),
		Canalization=genet.canalization(W, ..., mutsd=mutsd, replicates=replicates))
	return(ans)
}
