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

Rprof("prof.prof")
sim <- launchprogevol(myfile = param.file)
Rprof(NULL)

output.file <- ""
which.o <- which(cmd=="-o")
if (length(which.o > 0) && length(cmd) > which.o) {
	output.file <- cmd[which.o+1]
}

write.table(sim, file=output.file, sep="\t", quote=FALSE, row.names=FALSE)

