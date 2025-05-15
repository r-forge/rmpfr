### FIXME: Should become a short *wrapper* calling a more general
##  ~~~~~ plotRelErr <- function(x, FUN, fx = FUN(x), fxM = FUN(mpfr(x, precBits)), ...)
##        plotRelErr()  would become an exported  Rmpfr utility  function !
## <--> ./plot-factErr-def.R
##        ~~~~~~~~~~~~~~~~~~

##' Compute and Visualize (1--2 plots) the accuracy / relative error of gamma(xF-xG.R).
##' *Must* have the \pkg{Rmpfr} package attached, i.e., in `search()`.  { just for asNumeric() ?! }
##'
##' @title Compute and Visualize (In)accuracy of gamma(x)
##' @param x
##' @param gx numeric vector as x; defaults to `gamma(x)`.  Using `gamma(x)` computed
##'           alternatively, e.g., by an earlier R version allows interesing visual comparison
##' @param gamM an mpfr-vector as `x`; must be a correct `gamma(mpfr(x, precBits))`
##'             and is an argument only to allow skipping re-computation.
##' @param precBits integer number of bits to be for the high-accuracy \pkg{Rmpfr} gamma()
##' @param ch.Rversion the version of R used to computed `gx`
##' @param doFirst logical indicating if *two* plots (or just one) should be produced.
##' @param type
##' @param ylim1 y-limits of the first plot, i.e., for the (positive and negative) relative errors.
##' @param ylim2 y-limits of the 2n plot, i.e., for the *absolute relative errors.
##' @param cutBelow logical indicating if we should use a lower bound cutoff for the 2nd plot.
##' @param mar, mgp  graphical margin options, passed to `par()`.
##'
##' @return a numeric matrix, returned invisibly, with columns (x, gam, gamM, relErr).
##' @author Martin Maechler
p.gammEr <- function(x, gx = gamma(x), gamM = gamma(mpfr(x, precBits)),
                     precBits = 128, ch.Rversion = R.version.string,
                     doFirst = TRUE, type = "b",
                     vert = c(-10, 10, 50), # cutoffs in "old" gamma(): (-10, 10, 50)
                     ylim1=NULL, ylim2=NULL,
                     cutBelow = is.numeric(ylim2[1]),
                     mar = .1+c(3,4,4,4), mgp = c(1.6, .6, 0)) {
    rE <- asNumeric(sfsmisc::relErrV(gamM, gx))
    ##                       ^^^^^^^ (gx / gamM - 1) but dealing with Inf, NA, etc
    ## for plotExtr():
    fEps <-                          c(1/2, 1, 2,   10,   100,  1000)
    lEps <- eval(substitute(expression(E/2, E, 2*E, 10*E, 100*E,1000*E),
                            list(E = quote(epsilon[C]))))
    if(identical(ch.Rversion, R.version.string) && !capabilities("long.double") && !grepl("no-long", ch.Rversion))
        ch.Rversion <- paste0(R.version.string, "_no-long-double")
    plotExtr <- function() { # to be applied for both plots
        mtext(ch.Rversion, adj=1)
        abline(v = vert, col=3, lty=2)
        axis(1, at=vert, col=3, lty=2, col.axis=3)
        f <- fEps; if(!par("ylog")) { abline(h=0, lty=2);  f <- c(-rev(f),f) }
        f.eps <- f * .Machine$double.eps
        abline(h = f.eps, col="orange", lty=3)
        axis(4, at=f.eps[f.eps > 0], labels=lEps, col="orange", col.axis="orange",
             las = 2, cex=1/2, cex.axis = 3/4)
        legend("bottomleft", "x integer", col=2, text.col=2, pch=7, bty="n")
    }
    if(doFirst) { ## first plot : --------------------------
        op <- par(mfrow=2:1, mar=mar, mgp=mgp); on.exit(par(op))
        plot(x, rE, type=type, ylim = ylim1, cex=1/2,
             main="relative error of   gamma(x)")
        points(rE ~ x, subset = x %% 1 == 0, col=2, pch=7)
        plotExtr()
    } else { op <- par(mar=mar, mgp=mgp); on.exit(par(op)) }
    ## 2n plot : abs(<rel.errors>) --------------------------
    y <- abs(rE)
    if(cutBelow) {
        stopifnot(length(ym <- ylim2[1]) > 0, is.finite(ym))
        y[y < ym] <- 3/4*ym
    }
    plot(x, y, log = "y", type=type, ylim = ylim2, cex=1/2,
         yaxt = "n", # below use eaxis(2)
         ylab = quote(abs(rE)),
         main="|relative error| [log scale] of   gamma(x)")
    points(y ~ x, subset = x %% 1 == 0, col=2, pch=7)
    sfsmisc::eaxis(2)
    plotExtr()
    ## return matrix of results
    invisible(cbind(x, gam=gx, gamM=asNumeric(gamM), relErr = rE))
}
