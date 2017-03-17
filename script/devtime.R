#########################################################
# devtime.R
#
# Computes the "time to stability" of a W matrix
# 
# Copyright Arnaud Le Rouzic / CNRS 2016
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


script.dir <- normalizePath(dirname(parent.frame(2)$ofile))
if (is.null(script.dir))
	script.dir <- normalizePath(dirname(sys.frame(1)$ofile))

source(paste(script.dir, "netw.R", sep="/"))


devtime.M2 <- function(W, a=0.5, S0=rep(a, ncol(W)), thresh=0.01, max.steps=100, ...)
{
	stopifnot(thresh > 0, thresh < 1, max.steps > 0, ncol(W) > 0, ncol(W)==nrow(W))
	expression.time <- model.M2(W=W, a=a, S0=S0, steps=max.steps, full=TRUE, ...)$full
	diff.expression <- t(apply(expression.time, 1, diff))
	lessthanthresh <- abs(diff.expression) < thresh
	which(apply(lessthanthresh, 2, all))[1]
}
