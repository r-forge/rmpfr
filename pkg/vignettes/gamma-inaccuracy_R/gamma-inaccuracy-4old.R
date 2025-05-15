### This is to be run with R 4.2.1 (or a bit earlier) -- to produce pdf()s
### e.g.
"
R-4.2.1 CMD BATCH --no-save  gamma-inaccuracy-4old.R
"

require(Rmpfr)
require("sfsmisc")

## xF <- unique(sort(c(seq(-180, 180, by=1/16), outer(-180:180, c(-1/32, 1/64), "+"))))
## less nice binary numbers:
## xF <- unique(sort(c(seq(-177, 171, by=1/20), outer(-177:171, c(-.02, .01), "+"))))
sh <- function(x, iB=4L, iE=3L) { n <- length(x)
    cat(format(x[1:iB]),"...",format(x[(n-iE+1L):n]), "\n", sep=" ") }
## A mixture of nice (binary) k/16  and {nice decimal -- non-nice-binary} m/80 :
source("xF-xG.R")
## sh(i01, 8)
## stopifnot(all.equal(sixt, u01[iu16]))
sh(xG)
## Already sorted with "delta x" either 1/80 or 1/40:
## stopifnot(all.equal((1:2)/80, range(diff(xF))))
## stopifnot(all.equal((1:2)/80, range(diff(xG))))
length(xF) # 6803 (was 16575)

gammx <- gamma(xG)
facxF <- factorial(xF)
stopifnot(all.equal(gammx, facxF)) ## but not exactly:
all.equal(gammx, facxF, countEQ=TRUE, tol=0) ## Mean absolute diff: 2.062e-14

if(getRversion() <= "4.2.1") {
    gamRfile <- paste0("gamma_R-", as.character(getRversion()), ".rds")
    cat("Saving factorial(xF) to ", gamRfile,": ")
    saveRDS(structure(R.version = R.version.string,
                      cbind(xF = xF, xG = xG, fact.x = facxF, gamm.x = gammx)),
            file = gamRfile)
    cat("[Ok]\n")
}


if(interactive()) # not here typically
  withAutoprint({
    system.time(
      facMxF <- factorial(mpfr(xF, 256))
    ) # 0.757 sec
    relE <- asNumeric(facxF/facMxF - 1)
    cbind(xF, facxF, relE)
    (rx <- range(xF[abs(relE) < 1e-10])) # -171.5625  170.6125
    plot(relE ~ xF, type="l", lwd=1/2, ylim= c(-1,1)*2e-13, main="rel.error factorial(x)")
    cl <- "gray47"; abline(v=rx, lty=2, col=cl); axis(1, at=rx, col.axis=cl, col=cl)
    mtext(R.version.string, adj=1)
    ## zoom into x = [-9, 11]  <==>  x+1 in [-10, 10] { factorial(x) == gamma(x+1) }
    plot(relE ~ xF, subset=-11.5 < xF & xF < 9.5, type="l", lwd=1/2, ylim=c(-1,1)*2.5e-15)
    rx <- c(-11,9); cl <- adjustcolor(2,1/2)
    abline(v=rx, lwd=2, col=cl); axis(1, at=rx, col.axis=cl, col=cl)
    ##

    plot(pmax(1e-17, abs(relE)) ~ xF, subset=-11.5 < xF & xF < 9.5,
         type="l", log="y", ylim=c(1e-17, 8e-15))
    rx <- c(-11,9); cl <- adjustcolor(2,1/2)
    abline(v=rx, lwd=2, col=cl); axis(1, at=rx, col.axis=cl, col=cl)

})

## part of R CMD BATCH: proc.time()

