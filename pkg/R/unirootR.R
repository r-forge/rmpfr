### This is a translation of R_zeroin2 in ~/R/D/r-devel/R/src/appl/zeroin.c
### from  C to R  by John Nash,
### ---> file rootoned/R/zeroin.R of the new (2011-08-18) R-forge package rootoned
###
### Where John Nash calls it zeroin(), I call it unirootR()

##' Simple modification of uniroot()  which should work with mpfr-numbers
##' MM: uniroot() is in ~/R/D/r-devel/R/src/library/stats/R/nlm.R
##'
unirootR <- function(f, interval, ...,
		     lower = min(interval), upper = max(interval),
		     f.lower = f(lower, ...), f.upper = f(upper, ...),
                     extendInt = c("no", "yes", "downX", "upX"),
                     trace = 0,
		     verbose = as.logical(trace),
		     tol = .Machine$double.eps^0.25, maxiter = 1000,
                     check.conv = FALSE,
                     ## Rmpfr-only:
                     warn.no.convergence = !check.conv,
		     epsC = NULL)
{
    if(!missing(interval) && length(interval) != 2L)
	stop("'interval' must be a vector of length 2")
    ## For many "quick things", we will use as.numeric(.) but we do *NOT* assume that
    ## lower and upper are numeric!
    .N <- as.numeric
    if(lower >= upper) # (may be mpfr-numbers *outside* double.xmax)
	stop("lower < upper  is not fulfilled")
    if(is.na(.N(f.lower))) stop("f.lower = f(lower) is NA")
    if(is.na(.N(f.upper))) stop("f.upper = f(upper) is NA")

    ff <- f.lower * f.upper # used later also for mpfr-checking
    form <- if(inherits(ff, "mpfr")) function(x) format(x, drop0trailing=TRUE) else format
    Sig <- switch(match.arg(extendInt),
		  "yes" = NULL,
		  "downX"= -1,
		  "no"   =  0,
		  "upX"  =  1,
		  stop("invalid 'extendInt'; please report"))
    ## protect against later   0 * Inf  |--> NaN  and Inf * -Inf.
    truncate <- function(x) { ## NA are already excluded; deal with +/- Inf
        if(is.numeric(x))
            pmax.int(pmin(x, .Machine$double.xmax), -.Machine$double.xmax)
        else if(inherits(x, "mpfr") && is.infinite(x)) # use maximal/minimal mpfr-number instead:
            sign(x) * mpfr(2, .getPrec(x))^((1 - 2^-52)*.mpfr_erange("Emax"))
        else x
    }
    f.low. <- truncate(f.lower)
    f.upp. <- truncate(f.upper)
    doX <- (   is.null(Sig) && f.low. * f.upp. > 0 ||
	    is.numeric(Sig) && (Sig*f.low. > 0 || Sig*f.upp. < 0))
    if(doX) { ## extend the interval = [lower, upper]
	if(trace)
	    cat(sprintf("search in [%s,%s]%s", form(lower), form(upper),
			if(trace >= 2)"\n" else " ... "))
	Delta <- function(u) 0.01* pmax(1e-4, abs(u)) ## <-- FIXME? [= R's uniroot() for double]
        it <- 0L
	## Two cases:
	if(is.null(Sig)) {
	    ## case 1)	'Sig' unspecified --> extend (lower, upper) at the same time
	    delta <- Delta(c(lower,upper))
	    while(isTRUE(f.lower*f.upper > 0) &&
                  any(iF <- is.finite(c(lower,upper)))) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		if(iF[1]) {
		    ol <- lower; of <- f.lower
		    if(is.na(f.lower <- f(lower <- lower - delta[1], ...))) {
			lower <- ol; f.lower <- of; delta[1] <- delta[1]/4
		    }
		}
		if(iF[2]) {
		    ol <- upper; of <- f.upper
		    if(is.na(f.upper <- f(upper <- upper + delta[2], ...))) {
			upper <- ol; f.upper <- of; delta[2] <- delta[2]/4
		    }
		}
		if(trace >= 2)
		    cat(sprintf(" .. modified lower,upper: (%15g,%15g)\n",
				.N(lower), .N(upper)))
		delta <- 2 * delta
	    }
	} else {
	    ## case 2) 'Sig' specified --> typically change only *one* of lower, upper
	    ## make sure we have Sig*f(lower) <= 0 and Sig*f(upper) >= 0:
	    delta <- Delta(lower)
	    while(isTRUE(Sig*f.lower > 0)) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		f.lower <- f(lower <- lower - delta, ...)
		if(trace >= 2) cat(sprintf(" .. modified lower: %g\n", .N(lower)))
		delta <- 2 * delta
	    }
	    delta <- Delta(upper)
	    while(isTRUE(Sig*f.upper < 0)) {
		if((it <- it + 1L) > maxiter)
		    stop(gettextf("no sign change found in %d iterations", it-1),
			 domain=NA)
		f.upper <- f(upper <- upper + delta, ...)
		if(trace >= 2) cat(sprintf(" .. modified upper: %s\n", form(upper)))
		delta <- 2 * delta
	    }
	}
	if(trace && trace < 2)
            cat(sprintf("extended to [%s, %s] in %d steps\n", form(lower), form(upper), it))
    }
    if(!isTRUE(as.vector(sign(f.lower) * sign(f.upper) <= 0)))
	stop(if(doX)
	"did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0"
	     else sprintf("f() values at end points = (%s, %s) not of opposite sign",
                          form(f.lower), form(f.upper)))

    if(is.null(epsC) || is.na(epsC)) {
	## determine 'epsC'  ``the achievable Machine precision''
	## -- given the class of f.lower, f.upper
	if(is.double(ff)) epsC <- .Machine$double.eps
	else if(is(ff, "mpfr"))
	    epsC <- 2^-min(getPrec(f.lower), getPrec(f.upper))
	else { ## another number class -- try to see if getPrec() is defined..
	    ## if not, there's not much we can do
	    if(is(prec <- tryCatch(min(getPrec(f.lower), getPrec(f.upper)),
				   error = function(e)e),
		  "error")) {
		warning("no valid getPrec() for the number class(es) ",
			paste(unique(class(f.lower),class(f.upper)), collapse=", "),
			".\n Using double precision .Machine$double.eps.")
		epsC <- .Machine$double.eps
	    } else {
		epsC <- 2^-prec
		message("using epsC = %s ..", format(epsC))
	    }
	}
    }
    if(tol < epsC)
        warning(sprintf("tol (%g) < epsC (%g)  is rarely sensical,
 and the resulting precision is probably not better than epsC",
                        tol, epsC))

    ## Instead of the call to C code, now "do it in R" :
    ## val <- .Internal(zeroin2(function(arg) as.numeric(f(arg, ...)),
    ##			     lower, upper, f.lower, f.upper,
    ##			     tol, as.integer(maxiter)))
    a <- lower # interval[1]
    b <- upper # interval[2]
    fa <- f.lower # f(ax, ...)
    fb <- f.upper # f(bx, ...)
    if (verbose)
	cat(sprintf("Start zeroin: f(%g)= %g;  f(%g)= %g\n",
		    .N(a), .N(fa),  .N(b), .N(fb)))
    c <- a
    fc <- fa
    ## First test if we have found a root at an endpoint
    maxit <- maxiter + 2 # count evaluations as maxiter-maxit
    converged <- FALSE

    while(!converged && maxit > 0) { ##---- Main iteration loop ------------------------------

	if (verbose) cat("Top of iteration, maxit=",maxit,"\n")
	d.prev <- b-a
	## Distance from the last but one to the last approximation	*/
	##double tol.2;		       ## Actual tolerance		*/
	##double p;			## Interpolation step is calcu- */
	##double q;			## lated in the form p/q; divi-
	##				 * sion operations is delayed
	##				 * until the last moment	*/
	##double d.new;		## Step at this iteration	*/

	if(abs(fc) < abs(fb)) { ## Swap data for b to be the smaller
	    if (verbose) cat(sprintf("fc (=%s) smaller than fb\n", form(fa)))
	    a <- b
	    b <- c
	    c <- a	## best approximation
	    fa <- fb
	    fb <- fc
	    fc <- fa
	}
	tol.2 <- 2*epsC*abs(b) + tol/2
	if (verbose) cat("tol.2= ",.N(tol.2),"\n")
	d.new <- (c-b)/2 # bisection
	if (verbose) cat("d.new= ",.N(d.new),"\n")

	converged <- abs(d.new) <= tol.2  ||  fb == 0
	if(converged) {
	    if (verbose) cat("DONE! -- small d.new or fb=0\n")
	    ## Acceptable approx. is found :
	    val <- list(root=b, froot=fb, rtol = abs(c-b), maxit=maxiter-maxit)
	}
	else {
	    ## Decide if the interpolation can be tried	*/
	    if( (abs(d.prev) >= tol.2) ## If d.prev was large enough*/
	       && (abs(fa) > abs(fb)) ) { ## and was in true direction,
		## Interpolation may be tried	*/
		##    register double t1,cb,t2;
		if (verbose) cat("d.prev larger than tol.2 and fa bigger than fb --> ")

		cb <- c-b
		if (a == c) { ## If we have only two distinct points, linear interpolation
		    ## can only be applied
		    t1 <- fb/fa
		    p <- cb*t1
		    q <- 1 - t1
		    if (verbose) cat("a == c: ")
		}
		else { ## Quadric inverse interpolation*/
		    if (verbose) cat("a != c: ")
		    q <- fa/fc
		    t1 <- fb/fc
		    t2 <- fb/fa
		    p <- t2 * ( cb*q*(q-t1) - (b-a)*(t1-1) )
		    q <- (q-1) * (t1-1) * (t2-1)
		}
		if(p > 0) { ## p was calculated with the */
		    if (verbose) cat(" p > 0; ")
		    q <- -q ## opposite sign; make p positive */
		} else {    ## and assign possible minus to	*/
		    if (verbose) cat(" p <= 0; ")
		    p <- -p ## q				*/
		}
		if (p < 0.75*cb*q - abs(tol.2*q)/2 ## If b+p/q falls in [b,c]*/
		    && p < abs(d.prev*q/2)) {	## and isn't too large	*/
		    if (verbose) cat("p satisfies conditions for changing d.new\n")
		    d.new <- p/q ## it is accepted
		} else if(verbose) cat("\n")
		## If p/q is too large, then the
		## bisection procedure can reduce [b,c] range to more extent
	    }

	    if( abs(d.new) < tol.2) { ## Adjust the step to be not less than tolerance
		if (verbose) cat("d.new smaller than tol.2, adjusted to it.\n")
                d.new <- if(d.new > 0) tol.2 else -tol.2
	    }
	    a <- b
	    fa <- fb ## Save the previous approx. */
	    b <- b + d.new
	    fb <- f(b, ...)
            if (verbose) cat(sprintf("new f(b) = %s; ", form(fb)))
	    maxit <- maxit-1
	    ## Do step to a new approxim. */
	    if( ((fb > 0) && (fc > 0)) || ((fb < 0) && (fc < 0)) ) {
		if (verbose) cat(sprintf("make c to have sign opposite to b: f(c)=%s\n", form(fa)))
		## Adjust c for it to have a sign opposite to that of b */
		c <- a
		fc <- fa
            } else if(verbose) cat("\n")

	}## else not converged

    } ## end{ while(maxit > 0) } --------------------------------------------

    if(converged) {
	iter <- val[["maxit"]]
	if(!is.na(fb) &&  abs(fb) > 0.5*max(abs(f.lower), abs(f.upper)))# from John Nash:
	    warning("Final function magnitude seems large -- maybe converged to sign-changing 'pole' location?")
    } else { ## (!converged) : failed!
        if(check.conv)
            stop("no convergence in zero finding in ", iter, " iterations")
        ## else
	val <- list(root= b, rtol = abs(c-b))
	iter <- maxiter
	if(warn.no.convergence)
	    warning("_NOT_ converged in ", iter, " iterations")
    }
    list(root = val[["root"]], f.root = f(val[["root"]], ...),
	 iter = iter, estim.prec = .N(val[["rtol"]]), converged = converged)
}
