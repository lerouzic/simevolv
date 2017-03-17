#########################################################
# Mmatrix.R
#
# Computes and displays the mutation variance-covariance 
# matrix for a specific W matrix
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

Mmatrix.MonteCarlo.ij.run <- function(W, sdmut, i, j, what="mean", ...) {
	stopifnot(i > 0, i <= ncol(W), j > 0, j <= ncol(W))
	W[i,j] <- rnorm(1, mean=W[i,j], sd=sdmut)
	return(model.M2(W=W, ...)[[what]])
}

Mmatrix.MonteCarlo.ij <- function(W, sdmut=0.1, i, j, replicates=100, what="mean", ...) {
	dat <- do.call(rbind, lapply(1:replicates, function(x) Mmatrix.MonteCarlo.ij.run(W=W, sdmut=sdmut, i=i, j=j, what=what, ...)))
	ans <- list(ref=model.M2(W=W, ...)[[what]], data=dat, M=var(dat))
	class(ans) <- c("Mmatrix", class(ans))
	ans
}

Mmatrix.MonteCarlo <- function(W, sdmut=0.1, exclude.0=TRUE, replicates=100, what="mean", ...) {
	n <- ncol(W)
	ans <- list(ref=model.M2(W=W, ...)[[what]])
	for (i in 1:n) {
		for (j in 1:n) {
			if (W[i,j] != 0 || !exclude.0) {
				ans$data <- rbind(ans$data, Mmatrix.MonteCarlo.ij(W=W, sdmut=sdmut, i=i, j=j, replicates=replicates, what=what, ...)$data)
			}
		}
	}
	ans$M <- var(ans$data)
	class(ans) <- c("Mmatrix", class(ans))
	ans
}

plot.Mmatrix <- function(x, points=TRUE, matrix=TRUE, lim=NA, pch=1, lty=1, col.pch="gray", col.lty="red", level=0.95, ...) {
	n <- length(x$ref)
	stopifnot(n>1)
	mm <- matrix(0, ncol=n, nrow=n)
	mm[lower.tri(mm)] <- 1:(n*(n-1)/2)
	mm <- mm[-1,-n]
	layout(mm)
	
	par(mar=c(0,0,0,0)+0.1, oma=c(5,4,0,0))
	
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			xi <- x$data[,i]
			xj <- x$data[,j]
			plot(NULL, xlim=if(is.na(lim[1])) range(xi) else lim, ylim=if(is.na(lim[1])) range(xj) else lim, xlab="", ylab="", xaxt="n", yaxt="n")
			if (points) {
				points(xi, xj, col=col.pch, pch=pch)
			}
			if (matrix) {
				library(ellipse)
				lines(ellipse(x$M, level=level, which=c(i,j), centre=x$ref[c(i,j)]), lty=lty, col=col.lty)
				# plot eigenvectors? 
			}
			if (i == 1) { 
				axis(2)
				mtext(2, text=paste("Gene", j), padj=-3)
			}
			if (j == n) {
				axis(1)
				mtext(1, text=paste("Gene", i), padj=3)
			}
			points(x$ref[i], x$ref[j], pch=16, col="black", cex=2)			
		}
	}
	layout(1)
}
