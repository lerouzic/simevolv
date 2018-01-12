#########################################################
# canal.R
#
# runs canalization tests
#
# 
# Copyright Arnaud Le Rouzic / CNRS 2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

library(parallel)

library(R.utils, quiet=TRUE, warn=FALSE)
script.dir <- dirname(names(findSourceTraceback())[1])
source(paste(script.dir, "Wcanal.R", sep="/"))

options(mc.cores=detectCores()-2)

testcanal <- function(sizeW, muW=0, sdW=1, microSdE=0.1, mutsd=0.1, stab.thresh=1, density=1, target=rep(NA, sizeW), mindist=0.1*sizeW, max.tries=1e6, mat.replicates=100, can.replicates=100, fullD=FALSE, prog=FALSE, ...) {
    # if fullD, sizeW should be of size 1.
    stopifnot(!(length(sizeW) > 1 && fullD))
    mc.cores <- getOption("mc.cores", 2)
    
    pb <- if(prog && require(utils)) {txtProgressBar(min=1, max=mat.replicates, initial=1, style=3)} else NULL
    
    can.replicates <- round(can.replicates/mc.cores)*mc.cores # Optimal core use
    ll <- expand.grid(sizeW=sizeW, muW=muW, sdW=sdW, microSdE=microSdE, mutsd=mutsd)
    ll <- ll[rep(1:nrow(ll), mat.replicates), ]
    ans <- cbind(ll, do.call(rbind, mclapply(1:nrow(ll), function(ln) {
        if (ln %% mc.cores == 0 && !is.null(pb)) setTxtProgressBar(pb, ln)
        canW(W=makeRandW(ll$sizeW[ln], muW=ll$muW[ln], sdW=ll$sdW[ln], stab.thresh=stab.thresh, density=density, target=target, mindist=mindist, max.tries=max.tries, ...),
        replicates=can.replicates, 
        microSdE=ll$microSdE[ln],
        mutsd=ll$mutsd[ln], 
        fullD=fullD, ...)})))
    if (!is.null(pb)) { setTxtProgressBar(pb, mat.replicates);  close(pb)}
    return(ans)
} 

makeRandW <- function(n, muW=0, sdW=1, stab.thresh=1, density=1, target=rep(NA, n), mindist=0.1*n, max.tries=1e6, ...) {
    stopifnot(n > 0)
    stopifnot(density > 0 && density <= 1)
    stopifnot(length(target) == 1 || length(target) == n)
    stopifnot(sdW >= 0)
    stopifnot(stab.thresh > 0) 
    stopifnot(mindist > 0)
    stopifnot(max.tries > 0)
    
    num0 <- ceiling(n*n*(1-density))
    if (length(target) < n) target <- rep(target[1], n)
    
    bestW <- NA
    bestdist <- NA
    trynb <- 1
    repeat {
        currW <- matrix(rnorm(n*n, mean=muW, sd=sdW), ncol=n)
        currW[sample(1:(n*n), size=num0)] <- 0
        mm <- model.M2(W=currW, ...)
        if (mean(mm$var) < stab.thresh) {
            if (all(is.na(target))) {
                bestW <- currW    
            } else {
                di <- dist(rbind(target, mm$mean))
                if (is.na(bestdist) || di < bestdist) {
                    bestdist <- di
                    bestW <- currW
                }
            }
        }
        if (all(is.na(target)) && !all(is.na(bestW))) break
        if (!is.na(bestdist) && bestdist <= mindist) break
        if (trynb >= max.tries) {
            warning("makeRandW: max.tries (", max.tries, ") failed attempts.")
            break
        }
        trynb <- trynb + 1
    }
    return(bestW)
}

canW <- function(W, replicates=100, microSdE=0.1, mutsd=0.1, fullD=FALSE, index="avg", ...) {
    stopifnot(is.matrix(W))
    stopifnot(nrow(W) == ncol(W))
    if (any(is.na(W))) return(NA)
    n <- ncol(W)
    ans <- c(model.M2(W, ...)$mean, 
        dev.stability(W, index=index, ...), 
        homeostasis(W, index=index, microSdE=microSdE, ..., replicates=replicates),
        dev.canalization(W, index=index, ..., replicates=replicates),
        genet.canalization(W, index=index, mutsd=mutsd, ..., replicates=replicates),
        somatic.canalization(W, index=index, mutsd=mutsd, ..., replicates=replicates))
    nn <- c("Mean", "Stab", "Homeo", "EnvCan", "GenCan", "SomCan")
    names(ans) <- paste0(rep(nn, each=n), ".", 1:n)
    ans <- c(sapply(nn, function(ss) mean(ans[grepl(names(ans), pattern=ss)])), ans)
    names(ans)[1:length(nn)] <- nn
    if (fullD) ans else ans[1:length(nn)]
}
