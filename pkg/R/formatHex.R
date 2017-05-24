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

##' @title sprintf("%a", *)-like formatting of mpfr numbers
##' @param x mpfr-number vector
##' @param bits integer (scalar) specifing the desired number of bits ("binary digits")
##' @param style 1-character string specifying
##' @return character vector of same length as \code{x}
##' @author Martin Maechler
sprintfMpfr <- function(x, bits, style = "+") {
    stopifnot(length(style <- as.character(style)) == 1, nchar(style) == 1,
	      style %in% c("+", " "),
	      length(bits) == 1, bits %% 1 == 0)
    hexdigits <- 1L + (bits-1L) %/% 4L ## common to both branches
### TODO: For consistency, no longer use sprintf() for bits <= 52
### ----  currently "fails", e.g., in  mpfr(formatBin(mpfr(2, 60)))
    if(bits > 52) { # <== precBits > 53
	neg <- sign(x) == -1
	ff <- .mpfr2str(x, hexdigits + 1L, base = 16)  ## need +1
	isNum <- ff$finite	## ff$finite == is.finite(x)
	i0 <- ff$is.0	## == mpfrIs0(x)
	ex <- ff$exp ## the *decimal* exp : one too large *unless* x == 0
	r  <- ff$str # the mantissa, including "-" if negative
	Ex <- ex - 1L
	if(any(i0)) Ex[i0] <- ex[i0]
	if(!all(isNum)) ## "@Inf@", "@NaN@", ...
	    r[!isNum] <- gsub("@", '', r[!isNum], fixed=TRUE)
	if(any(i <- neg & isNum))
            ## r[i] <- sub("^-", "-0x", r[i])  wrongly gives e.g. "-0x.18"; want "-0x1.8"
	    r[i] <- paste0("-0x", substr(r[i], 2L, 2L), ".",
			   substring(r[i], 3L), "p")
	if(any(i <- !neg & isNum))
	    r[i] <- paste0(style, "0x", substr(r[i], 1L, 1L), ".",
			   substring(r[i], 2L), "p")
	r[isNum] <- paste0(r[isNum], c("", "+")[1+ (isNum & (Ex >= 0))], 4*Ex)
	r
    }
    else {
	nX <- as.character(hexdigits)
	sprintf(paste0("%", style, nX, ".", nX, "a"), x)
    }
}

formatHexInternal <- function(x, precBits = min(getPrec(x)), style = "+")
{
    if (is.numeric(x)) {
	precBits <- getPrec(x)
	x <- mpfr(x, precBits)
    }
    bindigits <- as.integer(precBits) - 1L
    hexdigits <- 1L + ((bindigits-1L) %/% 4L)
    ## hexdigits is the number of hex digits after the precision point
    structure(sprintfMpfr(x, bits=bindigits, style=style),
              ##---------
	      bindigits = bindigits,
	      hexdigits = hexdigits)
}# formatHexInternal

##___ ../man/formatHex.Rd ___
##           ~~~~~~~~~~~~
formatHex <- function(x, precBits = min(getPrec(x)), style = "+") {
    structure(formatHexInternal(x, precBits=precBits, style=style),
	      dim = dim(x), dimnames = dimnames(x),
              class = c("Hcharacter", "character"))
}

formatBin <- function(x, precBits=min(getPrec(x)), scientific = TRUE,
		      left.pad = "_", right.pad = left.pad, style = "+")
{
    H <- formatHexInternal(x, precBits=precBits, style=style)
    ## bindigits is number of binary digits after the precision point
    bindigits <- attr(H, "bindigits")
    hexdigits <- attr(H, "hexdigits")# *must* be correct = #{pure digits between "." and "p"}
    attributes(H) <- NULL
    S <- substr(H, 1, 1) # sign
    A <- substr(H, 4, 4)
    B <- substr(H, 6, 6+(hexdigits-1))
    ## assumes *always* an exponent "p" which is correct
    pow <- substr(H, 6+hexdigits+1, 1000000L)
    sB <- strsplit(B, "")
    rsB <- do.call(rbind, sB)
    hrsB <- HextoBin[rsB]
    dim(hrsB) <- dim(rsB)
    hrsBa <- apply(hrsB, 1, paste, collapse="")
    hrsBb <- substr(hrsBa, 1, bindigits)
    ## While this is a truncation,
    ## the mpfr conversion assures that
    ## only zero characters are truncated.
    if (!scientific) {
	powers <- as.integer(pow)
        Left <- -powers + max(powers, 2-precBits)
        Right <- powers - min(powers, precBits-1)
	D <- cbind(S, "0b", strrep(left.pad, Left),
		   A, hrsBb, strrep(right.pad, Right))
	D2 <- apply(D, 1, function(x) do.call(paste, list(x, collapse="")))
	ilft <- as.integer(max(Left) + min(powers)) + 4L
	res <- paste0(substr(D2,      1L,   ilft  ), ".",
		      substr(D2, ilft+1L, 1000000L))
    }
    else {
	res <- cbind(S, "0b", A, ".", hrsBb, "p", pow)
	res <- apply(res, 1, function(x) do.call(paste, list(x, collapse="")))
    }
    structure(res, dim = dim(x), dimnames = dimnames(x),
              class = c("Bcharacter", "character"))
}

print.Bcharacter <- function(x, ...) {
    print(unclass(x), quote=FALSE, right=TRUE, ...)
    invisible(x)
}
print.Hcharacter <- function(x, ...) {
    y <- unclass(x)
    ## FIXME? use  `attributes<-`(y, attributes(y)[...]) ?
    attr(y,"bindigits") <- NULL
    attr(y,"hexdigits") <- NULL
    print(y, quote=FALSE, right=TRUE, ...)
    invisible(x)
}


## MM Still thinks  that  formatMpfr() [ = format(<mpfr>) method !]
## --                      ----------
## should be made better _and_ be used *instead* of formatDec()


## RMH 2017-05-23, ~/R/MM/Pkg-ex/Rmpfr/formatDec-revised2.R :
formatDec <- function(x, precBits = min(getPrec(x)), digits=decdigits,
                      nsmall=NULL, scientific=FALSE,
                      style="+", decimalPointAlign = TRUE,
                      ...) {
    if (is.character(x)) x <- as.numeric(x)
    if (is.numeric(x)) x <- mpfr(x, precBits)
    else if (is.complex(x)) stop("complex 'x' are not supported in \"Rmpfr\" (yet)")
    decdigits <- ceiling(log(2^precBits, 10)) + 1
    ## NB: 'scientific' is ignored by format() here:
    chx <- format(x, digits=max(digits, decdigits), nsmall=nsmall,
                  scientific=scientific, style=style, ...)
    structure(if(decimalPointAlign) formatAlign(chx, ...) else chx,
	      dim = dim(x),
	      dimnames = dimnames(x),
	      class="noquote")
}

##' non exported utility  currently only used in  formatDec() :
formatAlign <- function(x, leftpad=" ", rightpad=leftpad, ...) {
  lr <- strsplit(x, ".", fixed=TRUE)
  l <- sapply(lr, `[`, 1) ## l left
  r <- sapply(lr, `[`, 2) ## r right
  r[is.na(r)] <- ""
  nl <- nchar(l)
  nr <- nchar(r)
  l.blank <- substr(paste0(rep( leftpad, max(nl)), collapse=""), 1L, max(nl) - nl)
  r.blank <- substr(paste0(rep(rightpad, max(nr)), collapse=""), 1L, max(nr) - nr)
  paste0(l.blank, l, ".", r, r.blank)
}



mpfr.Bcharacter <- function(x, precBits, scientific = NA, ...) {
    ## was scanBin()
    if (is.na(scientific)) ## we look for a "p" exponent..
        scientific <- any(grepl("p", x, fixed=TRUE))
    class(x) <- NULL
    noPrec <- (missing(precBits) || is.null(precBits))
    if (!scientific) {
        if (noPrec) {
	    ## Need to count the number of actual bits, and do it before converting "_" to "0".
	    ## Exclude  Inf, NaN, NA, ..
	    x. <- x[x %in% c("NaN", "NA", "Inf", "+Inf", "-Inf")]
	    ## "Unpad" {my first(!) example where a pipe is slightly more elegant}
	    x. <- gsub("_", "",
		       sub("_+$", "",
			   sub("^[-+]?0b", "",
			       sub(".", "", x, fixed=TRUE))), fixed=TRUE)
	    precBits <- max(nchar(x.), 1)
        }
        x <- gsub("_", "0", x) ## TODO: chartr(.......)
    }
    else if (noPrec) { ## scientific -- find precBits from looking at string:
	## assume a format such as "+0b1.10010011101p+1"
	no.p <- -1L == (i.p <- as.vector(regexpr("p", x, fixed=TRUE)))
	if(any(no.p)) ## no "p" - infer precision from others
	    i.p[no.p] <- round(mean(i.p[!no.p]))
	if(length(unique(i.p)) == 1L) ## all are the same
	    i.p <- i.p[1]
	precBits <- i.p - 5L
    }
    mpfr(x, base = 2, precBits=precBits, ...)
}

## A mpfr() method for "Hcharacter" .. but not quite -- called from mpfr() when appropriate
mpfr.Hcharacter <- function(x, precBits, ...) {
    class(x) <- NULL
    if (missing(precBits) || is.null(precBits)) {
	## assume a format such as "+0x1.745d17ap-4"
	no.p <- -1L == (i.p <- as.vector(regexpr("p", x, fixed=TRUE)))
	if(any(no.p))
	    i.p[no.p] <- round(mean(i.p[!no.p]))
	if( all(duplicated(i.p)[-1])) ## all are the same
	    i.p <- i.p[1]
	precBits <- 1 + (i.p - 6) * 4
    }
    mpfr(x, base = 16, precBits=precBits, ...)
}


`[.Bcharacter` <- `[.Hcharacter` <- ## == base :: `[.listof`
    function (x, ...) structure(NextMethod("["), class = class(x))

## Don't seem to get these to work correctly (at lost not easily):
## cbind.Bcharacter <- cbind.Hcharacter <-
##     function (...) structure(NextMethod("cbind"), class = class(..1))
## rbind.Bcharacter <- rbind.Hcharacter <-
##     function (...) structure(NextMethod("rbind"), class = class(..1))

