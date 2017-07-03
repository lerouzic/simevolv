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
