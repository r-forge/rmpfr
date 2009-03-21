#### All  coercion methods for the  "Rmpfr" classes


## was  toMpfr <-
mpfr <- function(x, precBits) {

    if(missing(precBits)) stop("must specify 'precBits'")
    stopifnot(length(precBits) == 1, precBits >= 2)
    ## libmpfr would exit (after good error message) for precBits == 1
    if(is.numeric(x))
        new("mpfr", .Call("d2mpfr1_list", x, precBits, PACKAGE="Rmpfr"))
    else if(is.character(x))
        stop("not yet implemented for character input")

    ## TODO: Make this work also with "character" x + optional argument 'base'


    else stop("invalid 'x'. Must be numeric or character")
}


setAs("numeric", "mpfr1", ## use default precision of 128 bits
      function(from) .Call("d2mpfr1", from, 128L, PACKAGE="Rmpfr"))
setAs("numeric", "mpfr", function(from) mpfr(from, 128L))
setAs("integer", "mpfr", function(from) mpfr(from, 128L))

setAs("mpfr", "numeric",
      function(from) .Call("mpfr2d", from, PACKAGE="Rmpfr"))
setMethod("as.numeric", "mpfr", function(x) .Call("mpfr2d", x, PACKAGE="Rmpfr"))
setMethod("as.integer", "mpfr", ## this is cheap; if needed do in libmpfr-code
	  function(x) as.integer(as.numeric(x)))

setAs("mpfr1", "numeric",  ## just for user-de-confusion :
      function(from) {
	  warning("coercing \"mpfr1\" via \"mpfr\" (inefficient)")
	  as(new("mpfr", list(from)), "numeric") })

setAs("mpfr1", "mpfr", function(from) new("mpfr", list(from)))


.mpfr2str <- function(x, digits = NULL) {
    stopifnot(is.null(digits) ||
	      (is.numeric(digits) && digits >= 0))
    ##	digits = NULL : use as many digits "as needed"
    .Call("mpfr2str", x, digits, PACKAGE="Rmpfr")
}

setMethod("format", "mpfr",
	  function(x, digits = NULL) {
	      stopifnot(is.null(digits) ||
			(is.numeric(digits) && digits >= 0))
	      ##  digits = NULL : use as many digits "as needed"

	      s.e <- .mpfr2str(x, digits)
	      ## The following could happen more efficiently
	      ## (but more error-prone) in C :
              ## TODO: still have 3 calls, 3 times R -> C:
              ## ----  .mpfr2str(), is.finite(),  mpfr.is.0()
              ## --> Return these results in *one* list !

	      ## (maybe) add decimal point
	      hasMinus <- sign(x) == -1
	      i. <- 1+hasMinus
              i0 <- mpfr.is.0(x)
	      if(!all(isNum <- is.finite(x))) {
		  r <- s.e$str
		  i. <- i.[isNum]
		  r[isNum] <- paste(substr(r[isNum], 1,i.),
				    substring(r[isNum], i.+1), sep = ".")
	      }
	      else
		  r <- paste(substr(s.e$str, 1,i.),
			     substring(s.e$str, i.+1), sep = ".")
	      ex <- s.e$exp # these are too large by one *unless* x == 0
	      if(any(i0))
		  ex[!i0] <- ex[!i0] - 1L
	      else ex <- ex - 1L
	      if(any(hasE <- isNum & ex != 0))
		  r[hasE] <- paste(r[hasE], as.character(ex[hasE]),
				   sep = "e")
	      r
	  })

setAs("mpfr", "character", function(from) format(from, digits=NULL))
