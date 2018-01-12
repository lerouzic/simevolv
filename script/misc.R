#########################################################
# misc.R
#
# Set of various helper functions.
# Nothing fancy here. 
# 
# Copyright Arnaud Le Rouzic / CNRS 2015
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

# mcsapply is a parallel-equivalent to sapply
"mcsapply" <- 
function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) 
{
    if (!require(parallel)) mclapply <- lapply
    FUN <- match.fun(FUN)
    answer <- mclapply(X = X, FUN = FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (!identical(simplify, FALSE) && length(answer)) 
        simplify2array(answer, higher = (simplify == "array"))
    else answer
}

# Implements the canalization index as defined in Rünneburger & Le Rouzic 2016 BMC Evol. Biol. 
# This is simply -log(V), where V is the variance in a disturbed phenotype.
"canIndex" <-
function (x, max.can=20)
{
    suppressWarnings(ans.tmp <- -log(x))
    ans.tmp[is.na(ans.tmp) | ans.tmp > max.can] <- max.can
    ifelse(is.na(x), NA, ans.tmp)
}

# Lightens a color
lighten <- function(color, factor=0.5){
    col <- col2rgb(color)
    col <- col+(255-col)*factor
    col <- rgb(t(col), maxColorValue=255)
    col
}


# redefine prcomp plots for customization
myplot.prcomp <- function(pca, PC1=1, PC2=2, pcol="gray", lcol="red", ...) {
    # Axes are oriented such as trait 1 has positive PCs 
    sig <- sign(pca$rotation[1,])
	x1 <- (t(pca$rotation)*sig)[PC1,]*8
    y1 <- (t(pca$rotation)*sig)[PC2,]*8
    pcax <- t(t(pca$x)*(sig))
	plot(pcax[,c(PC1,PC2)], xlab=paste0("PC",PC1), ylab=paste0("PC",PC2), pch=1, xlim=c(-1,1)+range(c(pcax[,PC1],x1)), ylim=c(-1,1)+range(c(pcax[,PC2],y1)), col=pcol, ...)
    arrows(0,0,x1=x1, y1=y1, col=lcol)
    text(x1, y1, label=rownames(pca$rotation), pos=ifelse(y1>0, 3, 1), col=lcol)
}

# principal component analysis for canalization
prcomp.analysis1 <- function(tab, what=c("Mean", "Stab", "Homeo", "EnvCan", "GenCan", "SomCan")) {
    layout(matrix(1:4, ncol=2))
    
    bygene.tab <- as.data.frame(lapply(what, function(ww) unlist(tab[,grepl(colnames(tab), pattern=paste0(ww, "."))])))
    colnames(bygene.tab) <- what
        
    myplot.prcomp(prcomp(tab[,what], scale=TRUE), main="Network average")
    myplot.prcomp(prcomp(tab[,what], scale=TRUE), main="Network average", PC2=3)    
    myplot.prcomp(prcomp(bygene.tab, scale=TRUE), main="By gene")
    myplot.prcomp(prcomp(bygene.tab, scale=TRUE), main="By gene", PC2=3)    
    layout(1)    
}

cor.analysis1 <- function(tab, what=c("Mean", "Stab", "Homeo", "EnvCan", "GenCan", "SomCan"), pcol="gray", lcol="black") {
    corplot <- function(x, y, capt=FALSE, captpos="bottom", wi=1, wj=1) {

        plot(x, y, col=pcol, xaxt="n", yaxt="n", xlab="", ylab="")
        mar <- 0.5*coef(lm(y ~ x))[2]+0.5*1/(coef(lm(x~y))[2])
        if (abs(cor(x, y)) > 0.2) abline(b=mar, a=mean(y)-mar*mean(x), col=lcol)
        text(min(x), max(y)-0.1*diff(range(y)), paste0("r=", round(cor(x, y), 2)), pos=4)
        if (capt) {
            if (captpos=="bottom") {
                axis(1)
                mtext(side=1, what[wi], line=5)
                axis(2)
            } else {
                axis(3)
                if (wi==1) mtext(side=3, what[wi], line=5)
                axis(4)
            }
        }
    }
    
    # The layout is a bit complex, but it helps to have it in the right order
    n <- length(what)
    mm <- matrix(NA, ncol=n, nrow=n)
    diag(mm) <- (n*n-n+1):(n*n)
    mm[lower.tri(mm)] <- 1:(n*(n-1)/2)
    mm <- t(mm)
    mm[lower.tri(mm)] <- (n*(n-1)/2+1):(n*(n-1))
    layout(mm)
    par(mar=rep(0.1, 4))
        
    bygene.tab <- as.data.frame(lapply(what, function(ww) unlist(tab[,grepl(colnames(tab), pattern=paste0(ww, "."))])))
    colnames(bygene.tab) <- what
    
    for (wi in 1:(n-1)) {
        for (wj in (wi+1):n) {
            corplot(tab[,what[wj]], tab[,what[wi]], capt=wj==wi+1, captpos="bottom", wj, wi)
        }
    }
    for (wi in 1:(n-1)) {
        for (wj in (wi+1):n) {
            corplot(bygene.tab[,what[wi]], bygene.tab[,what[wj]], capt=wj==wi+1, captpos="top", wi, wj)
        }
    }
    layout(1)    
}
