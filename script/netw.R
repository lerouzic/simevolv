#########################################################
# netW.R
#
# The version of the Wagner model used in 
# Rünneburger & Le Rouzic 2016 BMC Evol. Biol. 
#
# 
# Copyright Arnaud Le Rouzic / CNRS 2015-2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

sigma.M2 <- function(x, a) {
	1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
}

sigma.M2p <- function(x, lambda, mu) {
	1. / (1. + lambda * exp(-mu*x))
}

suppressMessages(library(compiler))
sigma.M2c <- cmpfun(sigma.M2p)

#~ library(Rcpp)

#~ cppFunction('
#~ NumericVector sigma_M2p(NumericVector x, double aam1, double l1am1) {
#~     NumericVector ans(x.size());
#~     // double aam1 = a*(1.-a);
#~     // double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

#~ cppFunction('
#~ NumericVector sigma_M2(NumericVector x, double a) {
#~     NumericVector ans(x.size());
#~     double aam1 = a*(1.-a);
#~     double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

model.M2 <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=20, measure=4, full=FALSE, varFUN=function(x) mean((x-mean(x))^2)) {
    lambda <- (1-a)/a
    mu <- 1/(a*(1-a))
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0 <- sigma.M2c((W %*% S0), lambda=lambda, mu=mu) 			
		sto[,i+1] <- S0
	}
	ans <- list()
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, varFUN)
	if (full) ans$full <- sto
	return(ans)
}

model.M2 <- cmpfun(model.M2)
