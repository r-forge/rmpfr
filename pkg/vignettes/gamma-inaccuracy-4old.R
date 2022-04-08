### This is to be run with R 4.1.3 -- to produce pdf()s
### e.g.
"
R-4.1.3 CMD BATCH --no-save  gamma-inaccuracy-4old.R
"

require(Rmpfr)
require("sfsmisc")

                    do.fig <- function() par(mar=c(5.1, 4.1, 1.1, 2.1))
##options(width = 75,
##	SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
## ...)

(Rver <- sub("[.]", "_", paste(R.version$major, R.version$minor, sep="_")))
## \SweaveOpts{engine=R,strip.white=true, width=8.5, height=6}
##                                             ^^^^        ^^
## Our custom graphics device:  vvvv        vv
pdfaCrop <- function(name, width=8.5, height=6, ...) {
    fn <- paste(paste0(name,"_R-", Rver), "pdf", sep = ".")
    if(FALSE)## debug
        cat("pdfaCrop: fn = ",fn,"; call:\n\t",deparse(match.call()),"\n")
    grDevices::pdf(fn, width = width, height = height, onefile=FALSE)# ...)
    do.fig() ## <<<<
    assign(".pdfaCrop.name", fn, envir = globalenv())
}
## This is used automagically :
pdfaCrop.off <- function() {
    dev.off()# for the pdf
    f <- get(".pdfaCrop.name", envir = globalenv())
    ## and now crop that file:
    pdfcrop <- "pdfcrop" # relying on PATH - fix if needed
    pdftex  <- "pdftex"  # relying on PATH - fix if needed
    system(paste(pdfcrop, "--pdftexcmd", pdftex, f, f, "1>/dev/null 2>&1"),
           intern=FALSE)
}



p.factEr <- function(x, precBits = 128, do2 = TRUE, type = "b", ylim1=NULL, ylim2=NULL) {
    fx <- factorial(x)
    facM <- factorial(mpfr(x, precBits))
    rE <- asNumeric(sfsmisc::relErrV(facM, fx))
    ##                       ^^^^^^^ (fx / facM - 1) but dealing with Inf, NA, etc

    fEps <-                          c(1/2, 1, 2,   10,   100)
    lEps <- eval(substitute(expression(E/2, E, 2*E, 10*E, 100*E),
                            list(E = quote(epsilon[C]))))
    plotExtr <- function() {
        mtext(R.version.string, adj=1)
        abline(v= -1 + c(10, 50), col=3, lty=2) # cutoffs in gamma() : (10, 50)
        f <- fEps; if(!par("ylog")) { abline(h=0, lty=2);  f <- c(-rev(f),f) }
        f.eps <- f * .Machine$double.eps
        abline(h = f.eps, col="orange", lty=3)
        axis(4, at=f.eps[f.eps > 0], labels=lEps, col="orange", col.axis="orange",
             las = 2, cex = 3/4)
        legend("bottomleft", "x integer", col=2, text.col=2, pch=7, bty="n")
    }
    if(do2) {
        op <- par(mfrow=2:1, mar=.1+c(2.5,3,3,2.5), mgp=c(1.6, .6, 0)); on.exit(par(op))
        plot(x, rE, type=type, ylim = ylim1, cex=1/2,
             main="relative error of   x! := factorial(x) := gamma(x+1)")
        points(rE ~ x, subset = x %% 1 == 0, col=2, pch=7)
        plotExtr()
    }

    plot(x, abs(rE), log = "y", type=type, ylim = ylim2, cex=1/2, yaxt="n",
         main="|relative error| [log scale] of   x! := factorial(x) := gamma(x+1)")
    points(abs(rE) ~ x, subset = x %% 1 == 0, col=2, pch=7)
    sfsmisc::eaxis(2)
    plotExtr()

    invisible(cbind(x, fac=fx, facM=asNumeric(facM), relErr = rE))
}

## plot-factEr
pdfaCrop("plot-factEr")
pfE1 <- p.factEr(seq(0, 49, by=1/8))
pdfaCrop.off()

## plot-factEr-full
pdfaCrop("plot-factEr-full")
pfE2 <- p.factEr(seq(0, 170, by=1/4))
pdfaCrop.off()

## plot-factEr-full-ylim
pdfaCrop("plot-factEr-full-ylim")
p.factEr(seq(0, 170, by=1/4), type="o",
         ylim1 = c(-1,1)*1e-14, ylim2 = c(.5e-17, 1e-13))
pdfaCrop.off()

## plot-factEr-neg-pos
pdfaCrop("plot-factEr-neg-pos")
pfEpm <- p.factEr(seq(-180, 180, by=1/8), type="l", ylim2= c(1e-17, 2e-13))
pdfaCrop.off()

(sI <- sessionInfo())

saveRDS(list(pfE1=pfE1, pfE2=pfE2, pfEpm=pfEpm, sInfo = sI),
        file = "gamma-inaccuracy_4old.rds")

proc.time()

