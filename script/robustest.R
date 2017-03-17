#########################################################
# robustest.R
#
# Computes and displays two robustness measurements
# used in Rünneburger & Le Rouzic 2016.
#
# * Indispensability test: to what extent a gene can be
#   extincted without influencing the rest of the newtork 
# 
# * Redundancy test: to what extent the expression of a specific
#   gene is influenced by the extinction of other network genes
# 
# Copyright Arnaud Le Rouzic / CNRS 2015
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################


sigma <- function(x, a) {
	1. / (1. + exp((-x/(a*(1.-a)))+log(1./a-1.)) )
}


model.genextinct <- function(W, extinct=NULL, a=0.2, S0=rep(a, nrow(W)), steps=100, measure=10) {
	# Runs the M2 model (cf netw.R) but forces a subset of genes to 
	# be extincted (0 expression). 
	# extinct can be NULL (equivalent to the full model), a scalar, or a vector. 
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0[extinct] <- 0
		S0 <- sigma(W %*% S0, a) 			
		sto[,i+1] <- S0
	}
	# Only returns the mean expression (no information about unstability)
	return (apply(sto[,(steps+1-measure):(steps+1)], 1, mean))
}

robustest <- function(W, a=0.2, S0=rep(a, nrow(W)), steps=100, measure=10) {
	ref <- model.genextinct(W=W, extinct=NULL, a=a, S0=S0, steps=steps, measure=measure)
	res <- sapply(1:nrow(W), model.genextinct, W=W, a=a, S0=S0, steps=steps, measure=measure)
	   # res is a square matrix. 
	   # columns: extincted genes
	   # lines: gene expressions
	ans <- list()
	ans$dist <- abs(t(res-ref))
		# Now the matrix is transposed.
		# ans$dist[i,j] is the effect of extincting gene i on the expression of gene j
	diag(ans$dist) <- NA
		# The interpretation of the diagonal is tricky (expression before forced extinction)
		# so better to remove it
	ans$infl <- apply(ans$dist, 1, mean, na.rm=TRUE)
	ans$robst <- apply(ans$dist, 2, mean, na.rm=TRUE)
		# These scores are called "Indispensability" and "Redundancy", respectively, 
		# in Rünneburger & Le Rouzic 2016. 
	return(ans)
}
