## NB: Assumes an R(-devel) version where  src/nmath/gamma.c has "extended" xmin
## --- when I wrote this, I had first,
##      # define xmin -178. // -177.56341258681965  was -170.5674972726612
## and then           -182.

has.Ldouble <- capabilities("long.double")
thisNam <- paste0("gamma-subnorm_",
                  if(has.Ldouble) "Ldouble" else "NO_long.double")

do.pdf <- FALSE
do.pdf <- !dev.interactive(TRUE)
if(do.pdf) { pFile <- paste0(thisNam, ".pdf"); pdf(pFile) }

x <- seq(-180.125, -176.5, length=4001)
## replace the negative integers by values very slightly to the left and right:
sum(ii <- x %% 1 == 0)# i.i := is integer
x <- sort(c(x[!ii], outer(x[ii], outer(c(-4:-1, 1:4), 2^-c(12,15)), `+`)))
length(x); cat(head(x), " .... ", tail(x, 2), "\n")
g <- gamma(x)
require(Rmpfr)
gM <- gamma(mpfr(x, 512))
## the "true", best re-normalized double precision values
B <- 1e20 ## use 1e20 for renormalization
gTe20 <- asNumeric(1e20 * gM)
pGam20 <- function(x,g,gT, ylim, ...) {
    plot (g  ~ x, type= "l", ylim=ylim, ylab = deparse(substitute(g)), ...)
    abline(h = 0, v = floor(min(x)) : ceiling(max(x)), lty=2, lwd=1/2, col="gray42")
    lines(gT ~ x, col = adjustcolor(2, 2/3), lwd=3, ...)
}
pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*5e-302, subset = x > -179)

##' compute fac * exp(lgamma(x)) *but* with correct sign *and* ensure no over/underflow from fac
gammaByLog <- function(x, fac = 1) {
    isNeg <- x < 0  &  (x %% 2) > 1 # gamma(x) < 0 for these
    sign <- 1 - 2*isNeg
    ## in gamma.c, have n := floor(x) - 1; sign <- if(n %% 2 == 0) -1 else 1
    sign * exp(log(fac) + lgamma(x))
}

cat("Has this R 'long double' ? ",  capabilities("long.double"), "\n")

## zoom in, now depending of double precision:
## =======
if(capabilities("long.double")) withAutoprint({

    pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*6e-303)
    pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*1e-303)
    ## ==> it seems to make sense to go all the way to -179.01
    lines(x, gammaByLog(x)*1e20, col = adjustcolor(4, 1/2), lwd=2, lty=2)# is *NOT* better
    all.equal(g, gammaByLog(x), tol=0) # see TRUE, just in case:
    table(g == gammaByLog(x)) # all TRUE

    ## This is much better: using log() *and* rescale to not get into subnormals:
    lines(x, gammaByLog(x, 1e20), col = adjustcolor(3, 2/3), lwd=3, lty=2)
    ## not perfect but still much better
    all.equal(gTe20, gammaByLog(x, 1e20), tol=0) # => 2.142e-12 ; then 6.3652e-14
    ## but then if we revert to "best" subnormals, there's no difference :
    all(g == gammaByLog(x, 1e20)/1e20)

}) else withAutoprint({ ## no long double .. ==> more severe rounding/underflow

    pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*40e-303)
    pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*20e-303, subset = x > -179.01)
    pGam20(x, g*1e20, gTe20, ylim = c(-1,1)*10e-303, subset = x > -179.01)
    lines (x, gammaByLog(x)*1e20,  col = adjustcolor(4, 1/2), lwd=2)# is better!!!
    ##===> we *should* use gammaByLog(x) here, instead of direct formula !!!!!!
    ## (in my version of R-devel gamma.c, we do since 2022-04-29, 17:51 !
    ##
    ## and this is even much better {but does not really help so much}:
    lines(x,  gammaByLog(x, 1e20), col = adjustcolor(3, 2/3), lwd=3, lty=2)
    all.equal(gammaByLog(x, 1e20)/1e20, g, countEQ=TRUE, tol=0) # 0.109.. !! -- now TRUE
    all.equal(gammaByLog(x),            g, countEQ=TRUE, tol=0) # (the same)

})

if(do.pdf) {
    dev.off()
    if(interactive())
        system(paste(getOption("pdfviewer"), pFile, "&"))
    else
        cat("finished plotting into file ", pFile, "\n")
}
