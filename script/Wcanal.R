#########################################################
# Wcanal.R
#
# Computes the canalization properties of a W matrix
# 
# Copyright Arnaud Le Rouzic / CNRS 2015-2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


library(R.utils)
script.dir <- dirname(names(findSourceTraceback())[1])
source(paste(script.dir, "netw.R", sep="/"))


########################### Quick documentation #######################
#
# ********** index ***********
#
# index = var: variance among mutants (or compared to the reference genotype, see 'reference')
# index = can: canalization score (-log(var))
# index = avg: average (absolute) effect of mutations
# index = pctx: proportion of replicates leading to an expression change larger than x%

canscale <- function(x) ifelse(is.na(x) | x < exp(-20), 20, -log(x))
pctscore <- function(index) {
    if (grepl(index, pattern="pct")) {
		pct <- as.numeric(substr(index, 4, nchar(index)))
		stopifnot(is.numeric(pct))
        return(pct/100)
    } else { return(NA) }
}
precond <- function(W, a=NULL, S0=NULL, index=NULL, which.i=NULL, which.j=NULL) {
    stopifnot(length(W) == 1 || is.matrix(W))
	stopifnot(nrow(W) == ncol(W))
    if (!is.null(a)) stopifnot(length(a) == 1 && (a > 0 || a < 1))
    if (!is.null(S0)) stopifnot(length(S0) == nrow(W))
    if (!is.null(S0)) stopifnot(all(S0 >= 0 & S0 <= 1))
    if (!is.null(index)) stopifnot(index %in% c("var","can","avg") || grepl(index, pattern="pct"))
 	if (!is.null(which.i)) stopifnot(all(which.i > 0 & which.i <= nrow(W)))
	if (!is.null(which.j)) stopifnot(all(which.j > 0 & which.j <= nrow(W)))
}

dev.stability <- function(W, index = "var", ...) {
	# Internal (genetic) stability at the end of the development
    
    precond(W=W, index=index)
    
    vF <- function(x) mean((x-mean(x))^2)
    if (index == "avg") 
        vF <- function(x) mean(abs(x-mean(x)))
    if (grepl(index, pattern="pct")) 
        vF <- function(x) mean(abs(x-mean(x)) >= pctscore(index))
	ans <- model.M2(W, varFUN=vF, ...)$var
    if (index == "can") {
        ans <- canscale(ans)
    }
    return(ans)
}

homeostasis <- function(W, index="var", a=0.5, S0=rep(a, nrow(W)), ..., microSdE=0.1, nb.affected.genes=ncol(W), delay=1, replicates=100) {
        # Influence of Gaussian microenvironmental perturbations at the end of the development
        # on the mean expression
        
    precond(W=W, a=a, S0=S0, index=index)
        
    ref <- model.M2(W, a=a, S0=S0, ...)$mean
    rr <- sapply(1:replicates, function(i) {
        newS0 <- ref
        affected.genes <- sample.int(length(newS0), nb.affected.genes)
        newS0[affected.genes] <- newS0[affected.genes] + rnorm(nb.affected.genes, 0, sd=microSdE)
        newS0[newS0 < 0] <- 0
        newS0[newS0 > 1] <- 1
        model.M2(W, a=a, S0=newS0, steps=delay, measure=1, ...)$mean })
	ans <- NULL
	rrt <- rr - ref
	if (index == "var" || index == "can") {
		ans <- apply(rrt, 1, function(x) mean(x^2))
		if (index == "can") ans <- canscale(ans)
	}
	if (index == "avg") {
		ans <- apply(abs(rrt), 1, mean)
	}
	if (grepl(index, pattern="pct")) {
		ans <- apply(abs(rrt) >= pctscore(index), 1, mean)
	}
	return(ans)
}

dev.canalization <- function(W, index="var", a=0.5, S0=rep(a, nrow(W)), ..., nb.affected.genes=ncol(W), replicates=100) {
        # Influence of the initial expression on the final mean expression

    precond(W=W, a=a, S0=S0, index=index)

    ref <- model.M2(W, a=a, S0=S0, ...)$mean
    rr <- sapply(1:replicates, function(i) {
        newS0 <- ref
        affected.genes <- sample.int(length(newS0), nb.affected.genes)
        newS0[affected.genes] <- runif(nb.affected.genes)
                model.M2(W=W, a=a, S0=newS0, ...)$mean
    })
	ans <- NULL
	rrt <- rr - ref
	if (index == "var" || index == "can") {
		ans <- apply(rrt, 1, function(x) mean(x^2))
		if (index == "can") ans <- canscale(ans)
	}
	if (index == "avg") {
		ans <- apply(abs(rrt), 1, mean)
	}
	if (grepl(index, pattern="pct")) {
		ans <- apply(abs(rrt) >= pctscore(index), 1, mean)
	}
	return(ans)
}

genet.canalization <- function(W, ..., index="var", which.i=1:nrow(W), which.j=1:ncol(W), mutsd=0.1, exclude.0=TRUE, replicates=100) {
	# This is the traditional canalization measurement (robustness to mutations)
	#
	# which.i and which.j specify the lines/columns that can be mutated
	# exclude.0 : disallow mutations of the 0s
	
    precond(W=W, index=index, which.i=which.i, which.j=which.j)
   
    n <- nrow(W)
	
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
		model.M2(W=myW, ...)$mean
	})
	ans <- NULL
	rrt <- rr - model.M2(W, ...)$mean

	if (index == "var" || index == "can") {
		ans <- apply(rrt, 1, function(x) sum(x^2)/length(x))
		if (index == "can") ans <- canscale(ans)
	}
	if (index == "avg") {
		ans <- apply(abs(rrt), 1, mean)
	}
	if (grepl(index, pattern="pct")) {
		pct <- as.numeric(substr(index, 4, nchar(index)))
		stopifnot(is.numeric(pct))
		ans <- apply(abs(rrt) >= pct/100, 1, mean)
	}
	return(ans)
}


somatic.canalization <- function(W, index="var", ..., which.i=1:nrow(W), which.j=1:ncol(W), mutsd=0.1, exclude.0=TRUE, delay=10, replicates=100) {
        # Influence of mutations at the end of the development
    
    precond(W=W, index=index, which.i=which.i, which.j=which.j)

	n <- nrow(W)
	
	mrow <- matrix(rep(1:n, n), nrow=n)
	mcol <- t(mrow)	
	avoided <- !(mrow%in%which.i) | !(mcol%in%which.j) | (W == 0 & exclude.0)
	stopifnot (sum(avoided) < length(W))   
    
    ref <- model.M2(W, ...)
        
	Wbox <- if (sum(!avoided) == 1) { rep(which(!avoided), replicates)}
			else {sample(which(!avoided), replicates, replace=TRUE)}
	Wdev <- rnorm(replicates, 0, sd=mutsd)    
    rr <- sapply(1:replicates, function(i) {
		myW <- W
		myW[Wbox[i]] <- W[Wbox[i]]+Wdev[i]
		model.M2(W=myW, S0=ref$mean, steps=delay, measure=5, ...)$mean
    })
	ans <- NULL
	rrt <- rr - ref$mean
	if (index == "var" || index == "can") {
		ans <- apply(rrt, 1, function(x) mean(x^2))
		if (index == "can") ans <- canscale(ans)
	}
	if (index == "avg") {
		ans <- apply(abs(rrt), 1, mean)
	}
	if (grepl(index, pattern="pct")) {
		ans <- apply(abs(rrt) >= pctscore(index), 1, mean)
	}
	return(ans)
}




####################### Additional functions ######################################

genet.canalization.ij <- function(W, index="var", ..., mutsd=0.1, replicates=100, what="mean", plot=FALSE, on.gene=1:nrow(W), zlim=c(-log(0.25), 20), text=TRUE) {
	# Computes the canalization score for every single regulatory element (i, j)
	dat <- expand.grid(1:nrow(W), 1:nrow(W))
	colnames(dat) <- c("i","j")
	dat <- cbind(dat, t(mapply(dat$i, dat$j, FUN=function(i, j) 
		canscale(genet.canalization(W=W, index=index, ..., which.i=i, which.j=j, 
			exclude.0=FALSE, mutsd=mutsd, replicates=replicates, 
			what=what)))
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
