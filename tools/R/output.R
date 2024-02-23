library(abind)


#Matrix extraction from output files
extract.P.mean <- function(tt, what="MPhen",  gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  ex
}

extract.matrix <- function(tt, what="CovPhen", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  if (sqrt(length(ex)) %% 1 != 0) stop("No way to make a square matrix out of ", length(ex), " elements.")
  matrix(ex, ncol=sqrt(length(ex)))
}

extract.fitness <- function(tt, what="Mfit", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  return(ex)
}

extract.opt <- function(tt, what="FitOpt1", gen=tt[nrow(tt), "Gen"]) {
  ex <- unlist(tt[tt[,"Gen"]==gen, grep(colnames(tt), pattern=what)])
  if (length(ex) == 0) stop("No match for ", what, " at generation ", gen,".")
  return(ex)
}

extract.P.matrix <- function(tt, gen=tt[nrow(tt),"Gen"]) {
  extract.matrix(tt, "CovPhen", gen=gen)
}

extract.M.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
  extract.matrix(tt, "MGenCov", gen=gen)
}

extract.W.matrix <- function(tt, gen=tt[nrow(tt), "Gen"]) {
  t(extract.matrix(tt, "MeanAll", gen=gen))
}

table.list.mean <- function(ll) {
	arr <- do.call(abind, c(ll, list(along=3)))
	rowMeans(arr, dims=2)
}

