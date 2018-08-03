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

### Contains functions to plot ToleranceCurves ("fitness" [or "distance"] f(environment) ) for GRN of the modified WM with plastic input.
	## Specifically, contains function to define(generate.choiceMatrix), generate(Get.dWij_TolCurve), and plot (Plot.dWij_TolCurve || Plot.dWij_Matrix) tolerance curves for "Knocked Out" GRNs.
	## The functions specified "dWij" are designed to reproduce deletions of a SINGLE INTERACTION, not of whole gene deletion (row and column)


# needed to assess different 'developments' (network kinetics) in different environments (values for 'input genes', aka 'plastic genes')
	# set the path to "plasticity_index" functions
source("./plasticity_index.R")



## Generates a boolean matrix indicating which interactions are to be set to 0 with "Get.dWij_TolCurve"

generate.choiceMatrix <- function(choiceMode, nRow, nCol=nRow, which.i, which.j, nodiag=TRUE, plastrate=c(1,rep(0, nRow-1)), plasticPromoterOverride=TRUE ) {
	# Generates a boolean matrix indicating which interactions to be set to 0 with "Get.dWij_TolCurve"
		# Dimensions of the matrix are provided in nRow and nCol. This function has not been designed to take into account non squares matrices.

	# argument "choiceMode" sets the way in which which.i/j are going to be interpreted.
		# "rowcol" will target the chosen columns of the whosen rows or vice-versa.
		# "target" will only delete the Wij with coordinates [which.i, which.j]. Coordinates should be vectors of same length which will be iterated upon
		# "all" needs no parametrization. It'll only take into account the "nodiag" argument and uses the "plastrate" argument (similar to "plasticity_index") to get the plastic gene(s)
			# note that it is similar to 'choiceMode="rowcol", which.i="all", which.j="all"'

## Argument tests
stopifnot(nRow==nCol) # Just in case : delete this line to use NON-square matrices... and see what happens !
if( !(choiceMode=="all" || choiceMode=="KO") ) {
	suppressWarnings(if(which.i=="all") which.i <- 1:nRow)
	suppressWarnings(if(which.j=="all") which.j <- 1:nCol)
	stopifnot( (is.numeric(which.i) && is.numeric(which.j)) )
}
stopifnot(sum(plastrate>0)>=1)	# La fonction d'obtention de la matrice de cibles devrait avoir du sens sans gène plastique. Quant à celle qui produit des courbes de tolerance...
	# si je teste un seul environnement sans gène plastique, ça peut avoir un sens... mais bon y'a ptetre plus simple pour faire des KOs.
if(!(choiceMode %in% c( "rowcol", "target", "all"))) { stop("Argument \"choiceMode\" must be one of \"rowcol\", \"target\", \"all\"") }

## Body
	choiceMatrix <- matrix(FALSE, nrow=nRow, ncol=nCol )

	if(choiceMode=="rowcol") {
		matRows <-  matrix(rep(1:nRow, nRow), nrow=nRow) # similar and faster to row(matrix(NA, nRow, nRow))
		targets <- (matRows %in% which.i) & (t(matRows) %in% which.j)
		#targets <- matrix(targets, nRow, nCol)		# ça marche tout aussi bien si j'ai pas besoin de definir une matrice de FALSES en amont... ai-je besoin de ce faire btw ?
		choiceMatrix[targets] <- TRUE

	} else if(choiceMode=="target") {
		stopifnot(length(which.i)==length(which.j))

		sapply(1:length(which.i), function(x) {
			choiceMatrix[which.i[x], which.j[x]] <- TRUE
			NewChoiceMatrix <- ( get("choiceMatrix", envir=parent.env(environment())) ) | ( get("choiceMatrix", envir=environment()) )
			assign("choiceMatrix", envir=parent.env(environment()), value=NewChoiceMatrix )
			} )

	} else if(choiceMode=="all") {
		choiceMatrix <- !choiceMatrix

	} else { stop("Argument \"choiceMode\" must be one of \"rowcol\", \"target\", \"all\"") }

	# Reset all diagonal elements to FALSE (usefull if they are already 0)
	if(nodiag) { diag(choiceMatrix) <- FALSE }
	# Set the row of plastic gene(s) to FALSE as they are not affected by reglation
	if(plasticPromoterOverride) {
		plastic <- which(plastrate>0)
		choiceMatrix[plastic,] <- FALSE
	}
return(choiceMatrix) }





# Splits an array to a list along a given dimension. Used in the following function ("Get.dWij_TolCurve")

split.along.dim <- function(array, dim)
	# splits an 'array' along any given 'dim'ension
	setNames(lapply(split(array, arrayInd(seq_along(array), dim(array))[, dim]),
		array, dim = dim(array)[-dim], dimnames(array)[-dim]),
		dimnames(array)[[dim]])





# Outputs an array or a list of 3 dimensions.
	# Dimensions 1 ("i") and 2 ("j") represent the position in which the tested GRN has been KO-ed.
	# 3rd dimension gives "index" (fitness or distance, refer to "plasticity_index") in all tested environment.
		# for the GRN knocked out for the appropriate interaction.

Get.dWij_TolCurve <- function(W, choiceMatrix, output="array", wildType=c(1,1), msg=FALSE, testedEnvironments=testedEnvs, plastrate=c(1,0,0,0,0,0), ...) {
	# This function will output the tolerance curve for matrix with the chosen interaction (Wij) deleted

	# Arguments in the ellipse are those of "plasticity_index"
	# "choiceMatrix" is a boolean matrix indicating which interactions to be knoked out (set to 0), typically produced with 'generate.choiceMatrix'
	# W is a GRN of the same dimensions of "choiceMatrix". It is supposed to possess a gene whose expression can be fixed and acts as an environmental signal.
		# this "plastic gene" is identified by the 'plastrate' argument passed to "plasticity_index".
		# while it isn't devoid of sense, this function will stop if no plastic gene is identified.

	# if (output="array")
		# "3 dimensional matrix". The first two dimensions are the Wij matrix (row=i, col=j).
		#		The third dimension is the environment. Use ans[i,j,] to get the tolerance curve of a given KO.
	# if (output="list")
		# Output is a list of list of vectors of numbers : [[i]] [[j]] [environment]
		#	[[i]][[j]] correspond to the corrdinates of the deleted interaction. [environment] contains the "index" ("fit" or "dist") in every tested environment ("signalEnv")

	# "wildType" should contain the coordinates of a position with no effect on development (diagonal or promoter of plastic gene, in default settings)
		# It can be set to NULL to avoid useless calculation if it is already obtained otherwise.
		# This is usefull to have the wild type development and compare it to knock Offs.

# Argument test :
stopifnot( output %in% c("array", "list") )
if(output=="list") { stopifnot( exists("split.along.dim") ) }
stopifnot( all(dim(W)==dim(choiceMatrix)) )
stopifnot( sum(plastrate > 0) >= 1)
stopifnot( is.numeric(testedEnvironments) && length(testedEnvironments)>=1 )

# Body :
if(!is.null(wildType) & all(is.numeric(wildType)) & length(wildType)==2) {
	choiceMatrix[wildType[1], wildType[2]] <- TRUE
	if(msg) { message(paste0("\tWild Type in coordinates [", wildType[1], ",", wildType[2], "]")) }
}
	Ans <- array(NA, c(dim(W), length(testedEnvironments)), dimnames=list( paste0("i", 1:nrow(W)) , paste0("j", 1:ncol(W)) , paste0("e", testedEnvironments) ) )

	KOcoord <- which(choiceMatrix, arr.ind=TRUE)
	RES <- apply(KOcoord, 1, function(x) {
			KoW <- W
			KoW[x[1],x[2]] <- 0
			plasticity_index(KoW, signalEnv=testedEnvironments, plastrate=plastrate, ...)
			})

stopifnot(nrow(KOcoord)==ncol(RES))
	for(mut in 1:ncol(RES)) {
		Ans[KOcoord[mut, "row"], KOcoord[mut, "col"], ] <- RES[,mut]
		}
stopifnot( length(which(!(is.na(Ans[,,1]))))==length(which(choiceMatrix))  && all( which(!(is.na(Ans[,,1])))==which(choiceMatrix) ) )

if(output=="array") { return(Ans)
} else if(output=="list") { 
	Ans <- split.along.dim(array=Ans, dim=3)
	return(Ans)
}}



# Arguments for plotmaking with the Plot.KO_TolCurve

plotArgs <- list(
		xlab="Environment",
		xlim=range(as.numeric(testedEnvs)),
		ylab="Fitness",
		ylim=c(0,1),
		main="Tolerance Curve",
		type="o"
		)

# Plot function. Conceived to plot curves from ko_env 3D arrays


Plot.dWij_TolCurve <- function(TolData, rowI, colJ, envE=testedEnvs, wildType=c(1,1), .legpos="bottomleft", .ncol=floor(curvesCount/10), ...,
			plotArgs=list(xlab="Environment",xlim=range(as.numeric(testedEnvs)),ylab="Fitness",ylim=c(0,1),main="Tolerance Curve",type="o") ) {
	# Plots data from a 3 dimensional array, typically obtained with Get.dWij_TolCurve

# Tests :
				## ça pourrait être des caractères ?
stopifnot( length(envE)>0 && is.numeric(envE) )
	suppressWarnings(if(rowI=="all") rowI <- 1:nrow(TolData))
	suppressWarnings(if(colJ=="all") colJ <- 1:ncol(TolData))
stopifnot( ( length(rowI)>0 && length(colJ)>0 ) && ( is.numeric(rowI) && is.numeric(colJ) ) )

stopifnot(is.array(TolData) && length(dim(TolData))==3)
	stopifnot( max(rowI) <= dim(TolData)[1] )
	stopifnot( max(colJ) <= dim(TolData)[2] )
	stopifnot( length(envE) <= dim(TolData)[3] )

if( sum( substring(dimnames(TolData)[[3]], 2) %in% as.character(envE) ) != length(envE) ) stop("Discrepancy between dimnames(TolData[[3]]) and envE")
# Except for verification purposes, envE is useless in this function.
	# It could be deduced from TolData ( substring(dimnames(TolData)[3][[1]], 2) )

stopifnot(is.list(plotArgs) )
# Body :
	do.call(plot, c(list("x"=NULL), plotArgs))

	curvesCount <<- 0
	WijCoord <<- c()

	for(rr in rowI) {
	for(cc in colJ) {
		curvesCount <<- curvesCount+1
		lines(x=as.numeric(envE), y=TolData[rr,cc,], col=rainbow(length(rowI)*length(colJ))[curvesCount], ...)
		WijCoord <<- rbind(WijCoord, c(rr,cc))
	}}
		if(!is.null(wildType)) {
		lines(x=as.numeric(envE), y=TolData[wildType[1],wildType[2],], col="black", pch=1, lwd=1.2, type="o")
		}
		if(curvesCount != nrow(WijCoord) ) { message("Some curves may be missing") }
		legend(.legpos, legend=paste0(WijCoord[,1], "-", WijCoord[,2] ), lty=1, col=rainbow(curvesCount), title="KO of W i-j", seg.len=0.5, ncol=.ncol)
}




plotArgs <- list(
		xlim=range(as.numeric(testedEnvs)),
		ylim=c(0,1.1),
		type="o",
		cex.axis=0.25,
		xlab=NA,
		ylab=NA,
		tck=-0.005,
		mgp=c(0,-0.001,0),	#c(3,1,0)=> default, position of c(Title, xticklab, yticklab)
		cex=0.3
		)

Plot.dWij_Matrix <- function(TolData, envE=testedEnvs, wildType=c(1,1), ...,
			plotArgs=list(xlim=range(as.numeric(testedEnvs)),ylim=c(0,1),type="l", 
				xlab=NA,ylab=NA,tck=-0.005,mgp=c(0,0,0),cex.axis=0.25, cex=0.3) ) {
		# Displays a matrix of plots showing the effect of each Wij deletion.

		# Arguments in the ellipse are those of "lines" (for the mutants)
			# PlotArgs ara passed to "plot" and "axis". As some arguments are meaningless to those functions

# Tests :
				## ça pourrait être des caractères ?
stopifnot( length(envE)>0 && is.numeric(envE) )
stopifnot(is.array(TolData) && length(dim(TolData))==3)
	stopifnot( length(envE) <= dim(TolData)[3] )

if( sum( substring(dimnames(TolData)[[3]], 2) %in% as.character(envE) ) != length(envE) ) stop("Discrepancy between dimnames(TolData[[3]]) and envE")
# Except for verification purposes, envE is useless in this function.
	# It could be deduced from TolData ( substring(dimnames(TolData)[3][[1]], 2) )

# This function is designed for square matrices.
stopifnot(ncol(TolData)==nrow(TolData))
roro <- nrow(TolData)

stopifnot(is.list(plotArgs) )

# Body :
	par(oma=c(0,0,1,0))
	layout(matrix(1:roro^2, roro, roro, byrow=FALSE))

	for( Wij in 1:roro^2) {
		thrdVec <- c() ; for(thrd in 0:(length(envE)-1)) { thrdVec <- append(thrdVec, Wij+400*thrd) }

		par(mar=c(0.25,0.25,0.25,0.25) )
		do.call(plot, c(list("x"=NULL), "bty"="n", "xaxt"="n", "yaxt"="n", plotArgs)) ; box(lwd=0.5, col="black")
#		if(plotArgs$"xaxt"=="n") { axis(1, at=as.numeric(envE), labels=TRUE, tck=plotArgs$"tck") }	# si je veux pouvoir controler de xaxt ou pas
	#changeY
		do.call(axis, c(list("side"=2), "at"=seq(plotArgs$"ylim"[1],plotArgs$"ylim"[2],0.1), list("labels"=FALSE), "lwd"=0.25, plotArgs))
		do.call(axis, list("side"=2, "at"=seq(plotArgs$"ylim"[1],plotArgs$"ylim"[2],0.1), "labels"=TRUE, "line"=0,
			"cex.axis"=plotArgs$"cex.axis","tcl"=plotArgs$"tck","mgp"=plotArgs$"mgp", "lty"=0) )
	#changeX
		do.call(axis, c(list("side"=1), "at"=as.numeric(envE), list("labels"=FALSE), "lwd"=0.25, plotArgs))
		do.call(axis, list("side"=1, "line"=-0.4, "labels"=TRUE, "at"=as.numeric(envE), "cex.axis"=plotArgs$"cex.axis", "tcl"=plotArgs$"tck", "mgp"=plotArgs$"mgp", "lty"=0))

		if( (anyNA(TolData[thrdVec])) ) {
			rect(xleft=plotArgs$"xlim"[1], plotArgs$"ylim"[1], plotArgs$"xlim"[2], plotArgs$"ylim"[2], col="darkgray")
		} else { 
		# WT
			lines(x=as.numeric(envE), y=TolData[wildType[1],wildType[2],], col="black", pch=1, lwd=1, type=plotArgs$"type", cex=plotArgs$"cex")
		# dWij
			lines(x=as.numeric(envE), y=TolData[thrdVec], ...)
		}
	}
}



testfunctions=FALSE
if(testfunctions) {
	tstMat <- oneMat.AdaptMANYind
	chcMat <- generate.choiceMatrix(choiceMode="all", nRow=20, which.i="all", which.j=c(1,3), nodiag=TRUE)

	TolData <- Get.dWij_TolCurve(W=tstMat, choiceMatrix=chcMat, output="array", wildType=c(1,1), msg=TRUE,
				testedEnvironments=testedEnvs, plastrate=c(1,0,0,0,0,0),
				opt=optcorr_Evol, index="fit", displayWarnings=FALSE, Ssel=fitSTR, Sprime=46000, outputPheno=FALSE)

		## taille minimale
	pdf(paste0("./testlayout.pdf"), width=24.5, height=24.5)
		Plot.dWij_Matrix(TolData=TolData, envE=testedEnvs, col="red", type="l", pch=4, plotArgs=plotArgs)

	dev.off()
}

