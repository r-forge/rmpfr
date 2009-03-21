require("Rmpfr")

## must take the *larger* of the two precisions:
stopifnot(format(mpfr(1, 60) / mpfr(7, 160)) ==
          "0.14285714285714285714285714285714285714285714285712")

(x <- mpfr(0:7, 100) / 7)
stopifnot( mpfr.is.0(x - x) ) # badly failed on 64-bit


eps2 <- 2 * .Machine$double.eps
eps8 <- 8 * .Machine$double.eps
stopifnot(all.equal(as.numeric(x+ 1L),
                    as.numeric(x)+1L, tol = eps2),
          (3 * x)/3 <= x,
          all.equal(as.numeric(x * 2L),
                    as.numeric(x + x), tol = 0))

## really want all.equal() *method* <<< FIXME
all.EQ <- function(x,y, tolerance = 2^-98, ...) {
    ## that this works, indeed, is proof that arithmetic  DOES  work
    environment(all.equal.numeric) <- environment()# so that our mean() is used!
    all.equal.numeric(x, y, ...)
}

u <- mpfr(0:17, 128)/17
two <- mpfr(2,100)
stopifnot(all.EQ(u ^ two, u ^ 2),
          identical(u ^ 2, u ^ 2L),
          all.EQ(two ^ u, 2 ^ u),
          identical(2 ^ u, 2L ^ u))

stopifnot(all.equal(as.numeric(x+1),
                    as.numeric(x)+1))

## When we compute with 100 bits,
## we should compare relative errors with  2^-100 :
prettyNum(format(abs((x+pi)-pi - x) / 2^-100), drop0 = TRUE)
