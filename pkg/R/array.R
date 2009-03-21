## From an "mpfr" object make an mpfr(Array|Matrix) :

setMethod("dim", "mpfrArray", function(x) x@Dim)
setMethod("dimnames", "mpfrArray", function(x) x@DimNames)

## 2 basic methods to construct "mpfr - arrays" ( mpfrArray | mpfrMatrix ) :

setMethod("dim<-", signature(x = "mpfr", value = "ANY"),
	  function(x, value) {
	      if(is.numeric(value) && all(value == (iv <- as.integer(value))))
		  new(if(length(iv) == 2) "mpfrMatrix" else "mpfrArray",
		      x, Dim = iv)
	  })

mpfrArray <- function(x, precBits, dim = length(x), dimnames = NULL) {
    new(if(length(dim) == 2) "mpfrMatrix" else "mpfrArray",
	.Call("d2mpfr1_list", x, precBits, PACKAGE="Rmpfr"),
	Dim = dim,
	DimNames = if(is.null(dimnames)) vector("list", length(dim))
		   else dimnames)
}

setAs("array", "mpfr", function(from)
      mpfrArray(from, 128L,
		dim = dim(from),
		dimnames = dimnames(from)))


setMethod("dimnames<-", signature(x = "mpfrArray", value = "ANY"),
	  function(x, value) {
	      if(!is.list(value)) stop("non-list RHS")
	      if(length(value) != length(x@Dim))
		  stop("RHS (new dimnames) differs in length from dim(.)")
	      x@DimNames <- value
	      x
	  })

setMethod("t", "mpfrMatrix",
	  function(x) {
	      d <- x@Dim; n <- d[1]; m <- d[2]
	      ## These are the indices to get the transpose of m {n x m} :
	      ## ind.t <- function(n,m)rep.int(1:n, rep(m,n)) + n*(0:(m-1))
	      x@.Data <- x@.Data[rep.int(1:n, rep(m,n)) + n*(0:(m-1))]
	      x@Dim <- c(m,n)
	      x@DimNames <- x@DimNames[2:1]
	      x
	  })
setMethod("t", "mpfr",
	  function(x) { # t(<n-vector>) |-->  {1 x n} matrix
	      r <- new("mpfrMatrix")
	      r@Dim <- c(1L, length(x))
	      r@.Data <- x@.Data
	      r
	  })

setMethod("as.vector", "mpfrArray", function(x) as(x, "mpfr"))
## a "vector" in  *one* sense at least ...

toNum <- function(from) {
    structure(.Call("mpfr2d", from, PACKAGE="Rmpfr"),
	      dim = dim(from),
	      dimnames = dimnames(from))
}

setAs("mpfrArray", "array", toNum)

setAs("mpfrMatrix", "matrix", toNum)


print.mpfrArray <- function(x, digits = NULL, ...) {
    stopifnot(is(x, "mpfrArray"), is.null(digits) || digits >= 2)
    ## digits = NA --> the inherent precision of x will be used
    n <- length(x)
    ch.prec <-
	if(n >= 1) {
	    rpr <- range(sapply(x, slot, "prec"))
	    paste("of precision ", rpr[1],
		   if(rpr[1] != rpr[2]) paste("..",rpr[2]), " bits")
	}
    cl <- class(x)
    p0 <- function(...) paste(..., sep="")
    cat(p0("'",cl,"'"), "of dim(.) = ",
        p0("(",paste(x@Dim, collapse=", "),")"),
        ch.prec, "\n")
    if(n >= 1) {
        ## FIXME: probably really need a 'format' method for these
        fx <- format(x, digits=digits)
        dim(fx) <- dim(x)
        dimnames(fx) <- dimnames(x)
	print(fx, ..., quote = FALSE)
    }
    invisible(x)
}
setMethod(show, "mpfrArray", function(object) print.mpfrArray(object))

## FIXME : should happen in C, where we could "cut & paste" much of
## -----  do_matprod() and matprod() from ~/R/D/r-devel/R/src/main/array.c
##/* "%*%" (op = 0), crossprod (op = 1) or tcrossprod (op = 2) */
.matmult.R <- function(x,y, op = 0)
{
    if(!(is.numeric(x) || is(x,"mpfr")))
        stop("'x' must be numeric of mpfr(Matrix)")
    sym <- missing(y)
    if (sym && (op > 0)) y <- x
    else if(!(is.numeric(y) || is(y,"mpfr")))
        stop("'y' must be numeric of mpfr(Matrix)")
    ldx <- length(dx <- dim(x))
    ldy <- length(dy <- dim(y))
    ## "copy, paste & modify" from  do_matprod():
    if (ldx != 2 && ldy != 2) {		#* x and y non-matrices */
	if (op == 0) {
	    nrx <- 1L; ncx <- length(x)
	} else {
	    nrx <- length(x); ncx <- 1L
	}
	nry <- length(y)
	ncy <- 1L
    }
    else if (ldx != 2) {		#* x not a matrix */
	nry <- dy[1]
	ncy <- dy[2]
	nrx <- ncx <- 0L
	if (op == 0) {
	    if (length(x) == nry) {	#* x as row vector */
		nrx <- 1L
		ncx <- nry # == length(x)
	    }
	    else if (nry == 1) {	#* x as col vector */
		nrx <- length(x)
		ncx <- 1L # == nry
	    }
	}
	else { #* crossprod */
	    if (length(x) == nry) {	#* x is a col vector */
		nrx <- nry # = length(x)
		ncx <- 1L
	    }
	}
    }
    else if (ldy != 2) {		#* y not a matrix */
	nrx <- dx[1]
	ncx <- dx[2]
	nry <- ncy <- 0L
	if (op == 0) {
	    if (length(y) == ncx) {	#* y as col vector */
		nry <- ncx # = length(y)
		ncy <- 1L
	    }
	    else if (ncx == 1) {	#* y as row vector */
		nry <- 1L # = ncx
		ncy <- length(y)
	    }
	}
	else {
	    if (length(y) == nrx) {	#* y is a col vector */
		nry <- nrx # = length(y)
		ncy <- 1L
	    }
	}
    }
    else {				#* x and y matrices */
	nrx <- dx[1]
	ncx <- dx[2]
	nry <- dy[1]
	ncy <- dy[2]
    }
    ##* nr[ow](.) and nc[ol](.) are now defined for x and y */

    z <- new("mpfrMatrix")
    z0 <- as(0, "mpfr")

    if (op == 0) { ## %*%
	if (ncx != nry) stop("non-conformable arguments")

        z@Dim <- c(nrx, ncy)
        z@.Data <- vector("list", nrx*ncy)
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
            for(i in 1:nrx)
                for (k in 0L:(ncy - 1L)) {
                    sum <- 0
                    for (j in 0L:(ncx - 1L))
                        sum <- sum + x[i + j * nrx] * y[1L+ j + k * nry]
                    z[i + k * nrx] <- sum
                }
        else        #/* zero-extent operations should return zeroes */
            for(i in seq_len(nrx*ncy)) z[i] <- z0
    }
    else if (op == 1) { ## crossprod() :  x' %*% y
	if (nrx != nry) stop("non-conformable arguments")

        z@Dim <- c(ncx, ncy)
        z@.Data <- vector("list", ncx*ncy)
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
            for(i in 0L:(ncx - 1L))
                for (k in 0L:(ncy - 1L)) {
                    sum <- 0
                    for (j in 1L:nrx)
                        sum <- sum + x[j + i * nrx] * y[j + k * nry]
                    z[1L +i + k * ncx] <- sum
                }
        else
            for(i in seq_len(ncx*ncy)) z[i] <- z0

    }
    else { ## op == 2 :  tcrossprod() :  x %*% y'
	if (ncx != ncy) stop("non-conformable arguments")

        z@Dim <- c(nrx, nry)
        z@.Data <- vector("list", nrx*nry)
        if (nrx > 0 && ncx > 0 && nry > 0 && ncy > 0)
            for(i in seq_len(nrx))
                for (k in 0L:(nry - 1L)) {
                    sum <- 0
                    for (j in 0L:(ncx - 1L))
                        sum <- sum + x[i + j * nrx] * y[1L +k + j * nry]
                    z[i + k * nrx] <- sum
                }
        else
            for(i in seq_len(nrx*nry)) z[i] <- z0

    }
    z
}

setMethod("%*%", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
          function(x,y) .matmult.R(x,y, op= 0))
## "FIXME"?  'ANY' is a bit much; may give uncomprehensible error messages
##       advantage: things will work with "Matrix" class matrices,
setMethod("%*%", signature(x = "mpfrMatrix", y = "array_or_vector"),
          function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "array_or_vector", y = "mpfrMatrix"),
          function(x,y) .matmult.R(x,y, op= 0))
## Even these (vectors!) -- if that's ok, drop the above two --
## FIXME --- ditto for  crossprod / tcrossprod :
setMethod("%*%", signature(x = "mpfr", y = "array_or_vector"),
          function(x,y) .matmult.R(x,y, op= 0))
setMethod("%*%", signature(x = "array_or_vector", y = "mpfr"),
          function(x,y) .matmult.R(x,y, op= 0))


setMethod("crossprod", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
          function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "mpfr", y = "array_or_vector"),
          function(x,y) .matmult.R(x,y, op= 1))
setMethod("crossprod", signature(x = "array_or_vector", y = "mpfr"),
          function(x,y) .matmult.R(x,y, op= 1))

setMethod("tcrossprod", signature(x = "mpfrMatrix", y = "mpfrMatrix"),
          function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "mpfr", y = "array_or_vector"),
          function(x,y) .matmult.R(x,y, op= 2))
setMethod("tcrossprod", signature(x = "array_or_vector", y = "mpfr"),
          function(x,y) .matmult.R(x,y, op= 2))


.mpfrA.subset <- function(x,i,j, ..., drop) {
    nA <- nargs()
    if(getOption("verbose"))
        message(sprintf("nargs() == %d  mpfrArray indexing ... ", nA))

    r <- x@.Data
    if(nA == 2) ## A[i]
        return(new("mpfr", r[i]))
    ## else: nA != 2 : nA > 2 -
    dim(r) <- (dx <- dim(x))
    dimnames(r) <- dimnames(x)
    r <- r[i,j, ..., drop=drop]
    if(drop & is.null(dim(r)))
        new("mpfr", r)
    else {
        x@Dim <- if(is.null(dr <- dim(r)))
            rep.int(1L, length(dx)) else dr
        x@DimNames <- if(is.null(dn <- dimnames(r)))
            vector("list", length(dx)) else dn
        attributes(r) <- NULL
        x@.Data <- r
        x
    }
}

## "["
setMethod("[", signature(x = "mpfrArray", i = "ANY", j = "ANY", drop = "ANY"),
          .mpfrA.subset)

## this signature needs a method here, or it triggers the one for "mpfr"
setMethod("[", signature(x = "mpfrArray", i = "ANY", j = "missing",
                         drop = "missing"),
          .mpfrA.subset)

## "[<-" :
setReplaceMethod("[", signature(x = "mpfrArray", i = "ANY", j = "ANY",
				value = "ANY"),
	function(x,i,j,value) {
            ## FIXME: should only trigger
            ##        for A[i,j] / A[i,] / A[,j] but not A[i]
            nA <- nargs()
            if(nA == 4) {
                r <- x@.Data
                dim(r) <- dim(x)
                dimnames(r) <- dimnames(x)
                r[i,j] <- as(value, "mpfr")@.Data
                attributes(r) <- NULL
                x@.Data <- r
            } else {
                stop(sprintf("nargs() == %d  mpfrArray[i,j] <- value ... ", nA))
            }
            x })
