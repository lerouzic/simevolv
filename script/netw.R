sigma.M2 <- function(x, a) {
	1. / (1. + exp((-x/(a*(1.-a)))+log(1./a-1.)) )
}


model.M2 <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=100, measure=10, full=FALSE, decay=1.0, microSdE=0) {
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0 <- S0*(1-decay)+decay*sigma.M2((S0 %*% W + rnorm(length(S0), 0, sd=microSdE)), a) 			
		sto[,i+1] <- S0
	}
	ans <- list()
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, function(x) mean((x-mean(x))^2))
	if (full) ans$full <- sto
	return(ans)
}

model.raw <-  function(W, S0=rep(1, nrow(W)), steps=100, measure=10, full=FALSE, decay=1.0) {
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0 <- S0*(1-decay)+decay*(S0 %*% W)
		sto[,i+1] <- S0
	}
	ans <- list()
	ans$mean <- apply(sto[,(steps+1-measure):(steps+1)], 1, mean)
	ans$var <- apply(sto[,(steps+1-measure):(steps+1)], 1, function(x) mean((x-mean(x))^2))
	if (full) ans$full <- sto
	return(ans)
}

plot.dyn <- function(W, y1=0, y2=0, steps=100, measure=10, model=model.M2, ylim=c(0,1), ...) {
	W[is.na(W)] <- c(y1, y2)
	mod <- model(W, steps=steps, measure=measure, full=TRUE, ...)
	plot(NULL, ylim=ylim, xlim=c(1, steps+1), xlab="Steps", ylab="Expression")
	for (i in 1:nrow(mod$full)) {
		lines(mod$full[i,], col=i)
	}
	abline(v=steps-measure, lty=2)
	legend("topleft", lty=1, col=1:nrow(mod$full), legend=as.character(1:nrow(mod$full)))
}

plot.Wmatrix <- function(W, ...) {
	plot(NULL, xlim=0:1, ylim=0:1, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", ...)
	arrows(x0=seq(0,1,length.out=nrow(W)+1), y0=0, y1=1, code=0)
	arrows(x0=0, x1=1, y0=seq(0,1,length.out=nrow(W)+1), code=0)
	shft <- 1/(2*nrow(W))
	unk <- 0
	for (y in 1:nrow(W)) {
		text(x=shft+2*shft*(0:(nrow(W)-1)), y=1-shft-2*shft*(y-1), ifelse(is.na(W[y,]), { unk <- unk+1; paste0("y",unk)}, W[y,]))
	}
}

plot.image <- function(data, variable, main=variable, xlab="y1", ylab="y2", zlim=c(0,1), col=gray(seq(0.9,0.1,length.out=101)), nlevels=4, ...) {
	mm <- matrix(data[,variable], nrow=length(unique(data$y1)))
	#mm <- t(mm)[nrow(mm):1,]
	image(x=unique(data$y1), y=unique(data$y2), z=mm, main=main, xlab=xlab, ylab=ylab, col=col, zlim=zlim, ...)
	if (var(c(mm), na.rm=TRUE) > 1e-10)
		contour(x=unique(data$y1), y=unique(data$y2), z=mm, add=TRUE, col="white", nlevels=nlevels)
}

test.bigene <- function(W, density=11, range=c(-5,5), a=0.5, steps=100, measure=10, accept.var=FALSE) {
	# W should contain one NA per line (the dimensions to be tested)
	
	alleles <- seq(from=range[1], to=range[2], length.out=density)
	ans <- expand.grid(list(y1=alleles, y2=alleles))
	ans$mean1 <- rep(NA, nrow(ans))
	ans$mean2 <- rep(NA, nrow(ans))
	ans$var1 <- rep(NA, nrow(ans))
	ans$var2 <- rep(NA, nrow(ans))	
	for (i in 1:nrow(ans)) {
		W.tmp <- W
		W.tmp[is.na(W.tmp)] <- unlist(ans[i, 1:2])
		getmod <- model.M2(W=W.tmp, a=a, steps=steps, measure=measure)
		ans[i,3:6] <- c(getmod$mean, getmod$var)
	}
	if (accept.var) {
		layout(cbind(c(1,2,4), c(0,3,5)))
		
		plot.Wmatrix(W)
		plot.image(ans, "mean1")
		plot.image(ans, "var1", zlim=c(0, 0.251))
		plot.image(ans, "mean2")
		plot.image(ans, "var2", zlim=c(0, 0.251))
	} else {
		layout(cbind(1:3))
		ans[ans$var1 > 1e-5, "mean1"] <- NA
		ans[ans$var2 > 1e-5, "mean2"] <- NA
		
		plot.Wmatrix(W)
		plot.image(ans, "mean1")
		plot.image(ans, "mean2")
	}
	layout(1)
	invisible(ans)
}
