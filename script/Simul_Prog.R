#!/usr/bin/env Rscript

#####################################################
#
# Simul_Prog.R
#
# R clone of the C++ Simul_Prog simulation program
# So far, many functions and parameters are missing
#
######################################################

cmd <- commandArgs(trailingOnly=FALSE)
my.path <- cmd[grep(cmd, pattern="--file=")]
my.path <- dirname(strsplit(my.path, split='=')[[1]][2])

source(paste(my.path, "simprogFun.R", sep="/"))

which.p <- which(cmd=="-p")
stopifnot(length(which.p)==1, length(cmd) > which.p)
param.file <- cmd[which.p+1]

#Rprof("prof.prof")
sim <- launchprogevol(myfile = param.file)
#Rprof(NULL)

output.file <- ""
which.o <- which(cmd=="-o")
if (length(which.o > 0) && length(cmd) > which.o) {
	output.file <- cmd[which.o+1]
}

write.table(sim[[1]], file=paste(date, output.file, ".table", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
library(rlist)
list.save(sim[[2]], file=paste(date, output.file, ".initialpop.rds", sep=""))
list.save(sim[[3]], file=paste(date, output.file, ".finalpop.rds", sep=""))
