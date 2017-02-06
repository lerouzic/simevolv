#########################################################
# misc.R
#
# Set of various helper functions.
# Don't take it seriously, nothing fancy here. 
# 
# Copyright Arnaud Le Rouzic / CNRS 2015
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

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

"canIndex" <-
function (x, max.can=20)
{
    suppresswarnings(ans.tmp <- -log(x))
    ans.tmp[is.na(ans.tmp) | ans.tmp > max.can] <- max.can
    ifelse(is.na(x), NA, ans.tmp)
}
