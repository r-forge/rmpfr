## From: Hans W Borchers <hwborchers@googlemail.com>
## To: Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: optimizeR for Rmpfr
## Date: Sun, 3 Jun 2012 16:58:12 +0200

## This is from Hans' pracma package,
## /usr/local/app/R/R_local/src/pracma/R/golden_ratio.R
## but there's also fibonacci search,  direct1d, ....
optimizeR <- function(f, lower, upper, ..., tol = 1e-20,
		      method = c("GoldenRatio"),
		      extraBits = 30,
		      precBits = -log2(tol) + extraBits, maxiter = 1000)
{
    stopifnot(length(lower) == 1, length(upper) == 1, lower <= upper)
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    a <- if(!is(lower,"mpfr")) mpfr(lower, precBits = precBits)
	 else if(.getPrec(lower) < precBits) roundMpfr(lower, precBits)
    b <- if(!is(upper,"mpfr")) mpfr(upper, precBits = precBits)
	 else if(.getPrec(upper) < precBits) roundMpfr(upper, precBits)
    method <- match.arg(method)
    ## if(method == "GoldenRatio") {
    switch(method, "GoldenRatio" = { ## golden ratio
	phi <- 1 - (sqrt(mpfr(5, precBits = precBits)) - 1)/2
	x <- c(a, a + phi*(b-a), b - phi*(b-a), b)
	y2 <- f(x[2])
	y3 <- f(x[3])
	n <- 0; convergence <- TRUE
	while ((d.x <- x[3] - x[2]) > tol) {
	    n <- n + 1
	    if (y3 > y2) {
		x[2:4] <- c(x[1]+phi*(x[3]-x[1]), x[2:3])
		y3 <- y2
		y2 <- f(x[2])
	    } else {
		x[1:3] <- c(x[2:3], x[4]-phi*(x[4]-x[2]))
		y2 <- y3
		y3 <- f(x[3])
	    }
	    if (n > maxiter) {
		warning(sprintf("not converged in %d iterations (d.x = %g)",
				maxiter, as.numeric(d.x)))
		convergence <- FALSE
		break
	    }
	}
	xm <- (x[2]+x[3])/2
	fxm <- if (abs(f. <- f(xm)) <= tol^2) 0. else f.
	list(min = xm, objective = fxm, iter = n,
	     convergence = convergence, estim.prec=abs(d.x))
    },
	   stop(sprintf("Method '%s' is not implemented (yet)", method)))
}
