#### Conversions   bigz  <-> mpfr   // also bigq <--> mpfr

if(packageVersion("gmp") < "0.5.8")## <-> ../NAMESPACE
    is.matrixZQ <- function(x) !is.null(attr(x, "nrow"))

## The following code is experimental, hence the "." :

### FIXME: we go via character.. which is not really efficient.
### Directly in C, we'd need both Rmpfr and gmp's  C code (!)
### TODO(?:  gmp should "export" its C++ API ( -> inst/include/*.hh )
### and we should add  'LinkingTo: gmp' to DESCRIPTION and
###  then use C++ with "C" { ...} for those parts
.bigz2mpfr <- function(x, precB = NULL) {
    stopifnot(inherits(x, "bigz"))
    ..bigz2mpfr(x, precB)
}
## Fast, no-checking (and not exported) version:
..bigz2mpfr <- function(x, precB = NULL)
    ## precB: 4 == log2(16) = log(base)
{
    b <- 16L
    cx <- .as.char.bigz(x, b)
    if(is.null(precB)) precB <- 4L*nchar(cx)
    if(is.matrixZQ(x))
	new("mpfrMatrix", .Call(str2mpfr1_list, cx, precB, b, "N"),
	    Dim = dim(x))# "bigz" has no dimnames
    else
	new("mpfr", .Call(str2mpfr1_list, cx, precB, b, "N"))
}
setAs("bigz", "mpfr", function(from) ..bigz2mpfr(from))


as.bigz.mpfr <-
.mpfr2bigz <- function(x, mod=NA) {
    if(is.null(mod)) mod <- NA_integer_
    stopifnot(is(x, "mpfr"),
	      is.na(mod) || (length(mod) == 1L && is.numeric(mod)))
    dx <- dim(x)
    cx <- format(trunc(x), drop0trailing=TRUE)
    dim(cx) <- dx ## needed?? {should *not* be, as in base R!}
    ..as.bigz(cx, mod)
}


## Fast, no-checking (and not exported) version:
..bigq2mpfr <- function(x, precB = NULL) {
    N <- numerator(x)
    D <- denominator(x)
    if(is.null(precB)) {
        eN <- frexpZ(N)$exp
        eD <- frexpZ(D)$exp
        precB <- pmax(128L, eN + eD + 1L) # precision of result
    }
    ..bigz2mpfr(N, precB) / ..bigz2mpfr(D, precB)
}

.bigq2mpfr <- function(x, precB = NULL) {
    stopifnot(inherits(x, "bigq"))
    ..bigq2mpfr(x, precB)
}
setAs("bigq", "mpfr", function(from) ..bigq2mpfr(from))
