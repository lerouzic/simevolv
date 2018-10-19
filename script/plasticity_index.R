#########################################################
# plasticity_index.R
#
# Computes an indicator of plasticity of a W matrix
# 
# Copyright Andreas Odorico / Paris-Sud University 2018
# <andreas.odorico@gmail.com>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################



## Outputs for "-W"
	#~ => mathematical mistake !!


#~	sigma.M2 <- function(x, a) {
#~		1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
#~	}

## Correct function
sigma.M2 <- function(x, a) {
	1. / (1. + ((1./a)-1)*exp(-x/(a*(1.-a))))
}



### Model with plasticity
model.plasti <- function(W, plastrate=c(0), plastsignal=c(0), a=0.2, S0=rep(a, nrow(W)), steps=24, measure=4, full=FALSE) {
	measure=measure-1	# so that "measure" and "steps" both 'start counting' at 1
		# note that S0 = sto[,1], hence the modification
	stopifnot(sum(plastrate>0)==length(plastsignal))
	stopifnot(length(S0)==nrow(W))
	stopifnot( (length(which(plastrate>1)) | length(which(plastrate<0)) | length(which(plastsignal<0)) | length(which(plastsignal>1)) )==0) #rmq... on peut juste ne pas inputter de la merde
	plastic=which(plastrate!=0)

	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
		#Re-evaluate S0 from standard S0 with plasticity
		S0[plastic] <- S0[plastic] + ((plastsignal[plastic]-S0[plastic])*plastrate[plastic])
	St <- S0
	sto[,1] <- S0
	for (i in 1:steps) {
		St <- sigma.M2(W %*% St, a)
		if(length(plastic)!=0){
			St[plastic] <- St[plastic] + ((plastsignal[plastic]-St[plastic])*plastrate[plastic])
		}
		sto[,i+1] <- St
	}
	ans <- list()
	if (measure == 0) {
		ans$val <- sto[,(steps-measure)]
	} else {
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, function(x) mean((x-mean(x))^2))
	}
	if (full) ans$full <- sto
	return(ans)
}



# Provides an index of plasticity given a W matrix and a expression of correlation between optima
	# example : opt <- "x ; x ; 1-x ; ((1/2)+(x/2)) ; ((1/2)-(x/2)) ; 0.5 ; 0.5 ; 0"

		## Added multiple indexes. Old versions had no various "index"es. To make them work, use " index="fit" ". Didn't put it as a default argument since fitness of the average individual should not be used as it differs from the average fitness in the population.

plasticity_index <- function(W, opt, index, signalEnv=seq(0,1,0.05), Ssel=c(0, 10, 10, 10, 10, 10), Sprime=46000, outputPheno=FALSE, displayWarnings=TRUE, model=model.plasti, plastrate=c(1,0,0), ...) {
			# multiple indexes :
			#	"fit" uses strength of selection, optima and expression to assess a fitness score in each "signalEnv" value
			#	"dist" outputs the distance to the optima


	### Multiple check-ups. No specific warnings.
stopifnot(ncol(W)==nrow(W))
N_loc <- ncol(W)
N_env <- length(signalEnv)
stopifnot( sum(index %in% c("dist", "fit"))==1 )

	### 1) make a matrix of [signal, optima]
if(!is.character(opt)){ stop("opt should be an expression as a function of \'x\', given as a string. The optima of every genes should be separated by a \";\" ")  }
OPTIMA <- unlist(strsplit(x=opt, split=";"))
# OPTIMA will fill with 0 to the size of N_loc
while(length(OPTIMA) < N_loc) { OPTIMA[length(OPTIMA)+1] <- 0 ; if(displayWarnings){message("\nlength(opt) < nrow(W) : auto-filling optimum with 0")} }

OPTIMA <- sapply(1:length(OPTIMA), function(oo) sapply(signalEnv, function(x) eval(parse(text=OPTIMA[oo])) ) )
	# outputs a matrix of [N_env, gene_optima]

	### 2) Development with every environmental signal
# plastrate autofills with 0.
while(length(plastrate) < N_loc) { plastrate[length(plastrate)+1] <- 0 ; if(displayWarnings){message("\nlength(plastrate) < nrow(W) : auto-filling plasticity rate with 0")} }

Pheno <- sapply(1:N_env, function(ss) model(W=W, plastsignal=OPTIMA[ss,], plastrate=plastrate, ...) )

Pheno <- list( t(sapply(1:N_env, function(x) unlist(Pheno[1,x]) )), t(sapply(1:N_env, function(x) unlist(Pheno[2,x]) )) )
names(Pheno)=c("mean", "var")
	# Pheno is a list of matrices which give : [[1]]=> mean, [[2]]=> var ; for every combination of [signal,gene]

	# Strength of stabilizing selection autofills unspecified values to the size of N_loc with 0
while(length(Ssel) < N_loc) { Ssel[length(Ssel)+1] <- 0 ; if(displayWarnings){message("\nlength(Ssel) < nrow(W) : auto-filling strength of selection on expression with 0")} } 

if(index=="fit") {
	# Strength of selection on developmental stability autofills unspecified values to the size of N_loc with itself (the first value]
	while(length(Sprime) < N_loc) { Sprime[length(Sprime)+1] <- Sprime[1] ; if(displayWarnings){message("\nlength(Sprime) < nrow(W) : auto-filling strength of selection on stability with the first value")} } 

	# Fitness calculation
	ansVEC <- sapply(1:N_env, function(ee) { prod(exp(-Ssel*(Pheno$mean[ee,] - OPTIMA[ee,])^2)) * prod(exp(-Sprime*Pheno$var[ee,])) } )
	names(ansVEC) <- signalEnv
		# returns a vector of fitness in every "N_env" condition

} else if (index=="dist") {

	if(sum(Ssel!=0)==0) stop("No selected gene : can't calculate distance to optimum")
	if(sum(OPTIMA[1,-1]!=0) != sum(Ssel!=0) && displayWarnings  ) {message("\nThere might be an inconsistency between opt[-1] and Ssel\n",
							"This may be due to selection for null gene expression (gene expression extinction) ; proceeding with distance calculation\n")}

	#ansVEC <- colMeans(sapply(1:N_env, function(ee) { abs(Pheno$mean[ee,] - OPTIMA[ee,]) } ))
	ansVEC <- sapply(1:N_env, function(ee) { abs(Pheno$mean[ee,(Ssel!=0)] - OPTIMA[ee,(Ssel!=0)]) } )

		if(is.vector(ansVEC) && !is.matrix(ansVEC) ) { ansVEC <- t(as.matrix(ansVEC)) }
	colnames(ansVEC) <- signalEnv
	ansVEC <- colMeans(ansVEC)
		# returns a vector of distance to optimum in every "N_env" condition
}

# Determination of Return
if(outputPheno){
	ans <- list(ansVEC, Pheno$mean, Pheno$var)
	names(ans)=c(index, "mean", "var")
}else{
	ans <- ansVEC
}
return(ans)
}



#~ Used to develop a character as a labile character.
	#~(i.e, starting the development in a next environment from the last value of the previous environment)

labile.plasti <- function(W, signalEnv, outputPheno=TRUE, plastrate=c(1, rep(0,19)), a=0.2, firstSzero=delayedAssign("firstSzero", rep(a, N_loc)), ...) {
	# in the ellipse are the arguments of "model.plasti"
	stopifnot(exists("model.plasti"))
	stopifnot(ncol(W)==nrow(W))
	N_loc <- ncol(W)
stopifnot(length(plastrate)==N_loc)

	eval(firstSzero)
	Szero <- firstSzero

	ans <- list(c(), c(), c(), c())
	names(ans) <- c("envs", "means", "vars", "devs")

	randEnv <- signalEnv[sample(order(signalEnv))]
for(EE in 1:length(signalEnv)) {
	res <- model.plasti(W=W, plastrate=plastrate, plastsignal=randEnv[EE], S0=Szero, a=a, full=outputPheno, ...)
	Szero <- res$mean

	if(outputPheno) {
		naames <- c()
		for(nnn in 0:(ncol(res$full)-1)) { naames <- append(naames, paste0("step", nnn, "_env", randEnv[EE])) }
			colnames(res$full) <- naames
			ans$devs <- cbind(ans$devs, res$full)
	
		ans$means <- cbind(ans$means, res$mean)
			colnames(ans$means)[ncol(ans$means)] <- paste0("env",randEnv[EE])
		ans$vars <-  cbind(ans$vars, res$var)
			colnames(ans$vars)[ncol(ans$vars)] <- paste0("env",randEnv[EE])

		ans$envs <- append(ans$envs, randEnv[EE])
	}
}
if(!outputPheno) { ans <- list(res$mean, res$var) ; names(ans) <- c("lastMeans", "lastVars") } #else { names(ans$envs) <- ans$envs }
return(ans)
}

