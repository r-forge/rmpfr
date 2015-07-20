## sprintf("%+13.13a", x) ## hex digits after the hex point = 13

## precBits: double precision = 53 = 1 + 13*4


## conversion from Hex digits to binary sequences of digits
HextoBin <- c(
 "0"="0000",
 "1"="0001",
 "2"="0010",
 "3"="0011",
 "4"="0100",
 "5"="0101",
 "6"="0110",
 "7"="0111",
 "8"="1000",
 "9"="1001",
 "A"="1010",
 "B"="1011",
 "C"="1100",
 "D"="1101",
 "E"="1110",
 "F"="1111",
 "a"="1010",
 "b"="1011",
 "c"="1100",
 "d"="1101",
 "e"="1110",
 "f"="1111")

if(FALSE) {
## the code isn't using either of these inverses.
BintoHex <- names( HextoBin[1:16])
names(BintoHex) <- HextoBin[1:16]

Bintohex <- tolower(BintoHex)
}


formatHexInternal <- function(x, precBits = min(getPrec(x)), style = "+") {
    if (precBits > 53) {
	## rmh warning
	warning("precBits reduced to 53 (sprintf does not currently support precBits > 53)")
	precBits <- 53
    }
    bindigits <- as.integer(precBits) - 1L
    hexdigits <- 1L + ((bindigits-1L) %/% 4L)
    ## hexdigits is the number of hex digits after the precision point
    ## rmh(2015-07-02)  if (missing(bindigits) && !missing(hexdigits)) bindigits <- 4*hexdigits
    format <- paste0("%", as.character(style),
		     as.character(hexdigits), ".", as.character(hexdigits), "a")
    structure(sprintf(format, x),
	      ##----- FIXME:  Fails for precisions higher than double
	      bindigits = bindigits,
	      hexdigits = hexdigits)
}# formatHexInternal


formatHex <- function(x, precBits = min(getPrec(x)), style = "+") {
    structure(formatHexInternal(x, precBits=precBits, style=style),
	      dim = dim(x), dimnames = dimnames(x), class = "Hcharacter")
}


## TODO : mpfr() should work with "formatBin" class
formatBin <- function(x, precBits=min(getPrec(x)), scientific = TRUE,
		      left.pad = "_", right.pad = left.pad, style = "+")
{
    H <- formatHexInternal(x, precBits=precBits, style=style)
    ## bindigits is number of binary digits after the precision point
    bindigits <- attr(H, "bindigits")
    hexdigits <- attr(H, "hexdigits")
    S <- substring(H, 1, 1)
    A <- substring(H, 4, 4)
    B <- substring(H, 6, 6+(hexdigits-1))
    pow <- substring(H, 6+hexdigits+1)
    sB <- strsplit(B, "")
    rsB <- do.call(rbind, sB)
    hrsB <- HextoBin[rsB]
    dim(hrsB) <- dim(rsB)
    hrsBa <- apply(hrsB, 1, paste, collapse="")
    hrsBb <- substring(hrsBa, 1, bindigits)
    ## While this is a truncation,
    ## the mpfr conversion assures that
    ## only zero characters are truncated.
    if (!scientific) {

	left.pads <- paste(rep(left.pad, 60), collapse="")
	left.pads <- substring(left.pads, 0, 0:60)

	right.pads <- paste(rep(right.pad, 60), collapse="")
	right.pads <- substring(right.pads, 0, 0:60)

	powers <- as.numeric(pow)
	Left <- -powers + max(powers)
	Right <- powers - min(powers)
	if (max(abs(powers)) > length(left.pads))
	    warning("Shifted binary out of bounds.", call.=FALSE)
	D <- cbind(S, "0b", left.pads[Left+1], A, hrsBb, right.pads[Right+1])
	D2 <- apply(D, 1, function(x) do.call(paste, list(x, collapse="")))
	res <- paste(substring(D2, 1, max(Left)+min(powers)+4), ".",
		     substring(D2, max(Left)+min(powers)+4+1), sep="")
    }
    else {
	res <- cbind(S, "0b", A, ".", hrsBb, "p", pow)
	res <- apply(res, 1, function(x) do.call(paste, list(x, collapse="")))
    }
    structure(res, dim = dim(x), dimnames = dimnames(x), class = "Bcharacter")
}

print.Hcharacter <-
print.Bcharacter <- function(x, ...) {
    print(unclass(x), quote=FALSE, right=TRUE, ...)
    invisible(x)
}


## formatDec1 <- function(x, ...) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10)+1)
##   result <- format(x, digits=decdigits, ...)
##   class(result) <- class(H)
##   result
## }

## this is the version I sent to Martin
## formatDec2 <- function(x, ..., digits, nsmall) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10))
##   Hx <- scan(text=H, quiet=TRUE)
##   absHx <- abs(Hx)
##   minabsHx <- min(absHx[absHx != 0]) ## actual 0 is wanted here
##   logminabsHx <- log(minabsHx %% 1, 10)
##   nsmall <- if (is.infinite(logminabsHx))
##               (decdigits-1) - log(min(absHx[absHx != 0]), 10)
##             else
##               ceiling(-logminabsHx)
##   result <- format(Hx, digits=decdigits, nsmall=nsmall, ...)
##   dim(result) <- dim(x)
##   dimnames(result) <- dimnames(x)
##   class(result) <- class(H)
##   result
## }

formatDec <- function(x, precBits = min(getPrec(x)), digits=decdigits,
                      nsmall=NULL, scientific=FALSE, style="+", ...) {
    H <- formatHexInternal(x, precBits=precBits, style=style)
    ## precBits <- attr(H, "bindigits") + 1
    decdigits <- ceiling(log(2^precBits, 10)) + 1
    Hx <- scan(text=H, quiet=TRUE) # <- FIXME ??? parse(.) or just  as.numeric() ?
    Hx.range <- range(abs(Hx[Hx != 0])) ## excluding hard zero
    if(!missing(digits))
	digits <- max(digits, decdigits)
    res <-
	if (scientific) {
	    sprintf(paste0("%", style, 1, ".", as.character(digits), "e"), Hx)
	}
	else {
	    if(is.null(nsmall)) nsmall <- -floor(log(Hx.range[1], 10))
	    if (nsmall <= 0) nsmall <- max(nsmall + digits, 0) ## effective only with scientific=FALSE
	    format(Hx, digits=digits, nsmall=nsmall, ...)
	}
    structure(res,
	      dim = dim(x),
	      dimnames = dimnames(x),
	      class = "noquote")
}

## Attempted sprintf.  doesn't align with "f", so can't use here.
## formatDec3 <- function(x, ..., displaydigits=decdigits, fmttype=c("e","f")) {
##   H <- formatHexInternal(x, ...)
##   precBits <- attr(H, "bindigits") + 1
##   decdigits <- ceiling(log(2^precBits, 10))
##   Hx <- scan(text=H, quiet=TRUE)
##   absHx <- abs(Hx)

##   Hx.range <- range(abs(Hx[Hx != 0]))
##   fieldwidth <- ceiling(diff(log(Hx.range, 10)))
##   fmttype <- match.arg(fmttype)
##   precision <- log(Hx.range[1], 10)
##   precision <- if (precision >= 0) 0 else -floor(precision)
##   digits <- max(decdigits, displaydigits)
##   switch(fmttype,
##          "e"=format <- paste("%+", as.character(digits), ".", as.character(digits), "e", sep=""),
##          "f"=format <- paste("%+", as.character(digits), ".", as.character(precision), "f", sep="")
##          )

##   result <- sprintf(format, Hx)

##   dim(result) <- dim(x)
##   dimnames(result) <- dimnames(x)
##   class(result) <- class(H)
##   result
## }


### FIXME:  mpfr() should be smart enough for these:

## mpfr(<"Bcharacter">) should do this
## a mpfr() method for "Bcharacter" .. but not quite -- called from mpfr() when appropriate
mpfrBchar <- function(x, precBits, scientific = TRUE, ...) {
  ## was scanBin()
  class(x) <- NULL
  if (!scientific) {
    x <- gsub("_", "0", x) ## TODO: chartr(.......)
    if (missing(precBits)) {
      ## mm: why warning?  Rather just make this the explicit default ?
      ## rmh: no.  we need to count the number of actual bits, and do it before converting "_" to "0".
      precBits <- mpfr_default_prec()
      warning("Default precBits = ", precBits)
    }
  }
  if (missing(precBits) || is.null(precBits)) {
    x1 <- strsplit(x[1],"")[[1]]
    plocation <- charmatch("p", x1) # FIXME: fails if there is no "p"
    precBits <- nchar(paste0(substring(x[1], 4, 4),
			     substring(x[1], 6, plocation-1)))
  }
  mpfr(x, base=2, precBits=precBits, ...)
}

## mpfr(<"Hcharacter">) should do this
## a mpfr() method for "Hcharacter" .. but not quite -- called from mpfr() when appropriate
mpfrHchar <- function(x, precBits, ...) {
  class(x) <- NULL
  if (missing(precBits) || is.null(precBits)) {
    x1 <- strsplit(x[1],"")[[1]]
    plocation <- charmatch("p", x1) # FIXME: fails if there is no "p"
    precBits <- 1 + nchar(substring(x[1], 6, plocation-1))*4
  }
  mpfr(x, base=16, precBits=precBits, ...)
}
