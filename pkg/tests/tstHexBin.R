require(Rmpfr)
sessionInfo()

nums <- c(0, .625, .1, .3, .3125, .5, .6, (.3+.6), .9, 1, 3.3125)
nums9 <- mpfr(nums, precBits = 9)

n5.b2 <- mpfr("101", base=2)
stopifnot(identical(n5.b2, mpfr(5, precBits=3)),
          identical(n5.b2, mpfr("0b101", base=2)))
if(FALSE)## FIXME (?? : 0b101  could also be a Hex number with no initial '0x')
identical(n5.b2, mpfr("0b101"))
## we *could* say that  anything starting with '0b' is binary, anything with '0x' is hexadecimal


### mpfrBchar (was 'scanBin') :

##' Check the inverse of formatBin(), i.e., mpfr() working correctly
chkInv.fBin <- function(x, ...) {
    stopifnot(is(x, "mpfr"))
    nb <- formatBin(x, ...)
    xx <- mpfr(nb)
    ## Inverse should work {apart from 0: that is not uniquely represented in MPFR!}
    stopifnot(identical(mpfr.is.0(x ) -> i0,
                        mpfr.is.0(xx)),
              identical(x[!i0], xx[!i0]))
    invisible(nb)
}
(nums9bc  <- chkInv.fBin(nums9))
(nums9bF <- chkInv.fBin(nums9, scientific=FALSE)) # "F": Fixed format (as in Fortran)

## higher precision, 0,Inf, sign change:
(i12 <- 1 / mpfr(c(-2:12, Inf), 64))
(i12.50 <- roundMpfr(i12, precBits=50)) # "same", with 50 (< 53) bits
try({ ## FIXME !!
(nI.12    <- chkInv.fBin(i12   ))
(nI.12.50 <- chkInv.fBin(i12.50))
})

n9. <- Rmpfr:::mpfrBchar(nums9bc)
n9_ <- mpfr             (nums9bc)
## Inverse worked {apart from 0: it is not at all uniquely represented in MPFR!}
stopifnot(identical(mpfr.is.0(n9.), mpfr.is.0(n9_)),
          identical(nums9[-1], n9.[-1]))

mpfr(nums9bc, precBits=5)

(n9.5 <- mpfr(nums9bF, scientific=FALSE, precBits=5))
stopifnot(all.equal(n9.5, mpfr(nums9bF, precBits=5), tol=0))
stopifnot(all.equal(n9., n9.5, tol = 0.02), getPrec(n9.5) == 5)
mpfr(nums9bF, scientific=FALSE)
mpfr(nums9bF)

### mpfr() -> mpfrHchar (was 'scanHex') :
nums9hc <- formatHex(nums9)
stopifnot(all.equal(n9.,  mpfr(nums9hc),             tol = 0),
          all.equal(n9.5, mpfr(nums9hc, precBits=5), tol = 0))

### Rest from ../R/formatHex.R -- heavily modified

(m32 <- matrix(0:31, 8, 4, dimnames = list(0:7, c(0,8,16,24))))
FourBits <- mpfr(m32, precBits=4)
FiveBits <- mpfr(m32, precBits=5)
stopifnot(all.equal(m32, asNumeric(FiveBits), tol=0),
          all.equal(m32, asNumeric(FourBits), tol= 0.05))

## HH:::formatHexInternal(FourBits)
formatHex(FourBits)
formatHex(FourBits, precBits=5)
formatHex(FourBits, precBits=6)
formatBin(FourBits)
formatBin(FourBits, precBits=5)
formatBin(FourBits, scientific=FALSE)
formatBin(FourBits, scientific=FALSE, precBits=5)
formatDec(FourBits)
formatDec(FourBits, precBits=5)
formatDec(FourBits, precBits=7)

## HH:::formatHexInternal(FiveBits)
formatHex(FiveBits)
formatHex(FiveBits, precBits=5)
formatHex(FiveBits, precBits=6)
formatBin(FiveBits)
formatBin(FiveBits, precBits=5)
formatBin(FiveBits, precBits=6)
formatBin(FiveBits, scientific=FALSE)
formatBin(FiveBits, scientific=FALSE, precBits=6)
formatDec(FiveBits)
formatDec(FiveBits, precBits=5)
formatDec(FiveBits, precBits=7)

TenPowers <- mpfr(matrix(10^(-3:5)), precBits=53)
row.names(TenPowers) <- -3:5
TenPowers
## HH:::formatHexInternal(TenPowers)
formatHex(TenPowers)
formatBin(TenPowers)
formatBin(TenPowers, scientific=FALSE)
formatDec(TenPowers)
formatDec(TenPowers, scientific=FALSE)
formatDec(TenPowers, precBits=54)
formatDec(TenPowers, precBits=54, scientific=FALSE)

TenFrac <- matrix((1:10)/10, dimnames=list(1:10, expression(1/x)))
TenFrac
formatHex(TenFrac)
formatBin(TenFrac)
formatBin(TenFrac, scientific=FALSE)
formatDec(TenFrac)

TenFrac9 <- mpfr(TenFrac, precBits=9)
TenFrac9
formatHex(TenFrac9)
formatBin(TenFrac9)
formatBin(TenFrac9, scientific=FALSE)
formatDec(TenFrac9)

stopifnot(getPrec(TenFrac ) == 53,
          getPrec(TenFrac9) ==  9)


Ten <- matrix(1:10 + 0.0, dimnames=list(1:10, "x")) ## + 0.0 forces double precision
Ten
formatHex(Ten)
formatHex(Ten, precBits=4)
formatBin(Ten)
formatBin(Ten, precBits=4)
formatBin(Ten, scientific=FALSE)
formatBin(Ten, scientific=FALSE, precBits=4)
formatDec(Ten)
formatDec(Ten, displaydigits=4)
formatDec(Ten, displaydigits=4, scientific=TRUE)

Ten4 <- mpfr(Ten, precBits=4)
Ten4
formatHex(Ten4)
formatBin(Ten4)
formatBin(Ten4, scientific=FALSE)
formatDec(Ten4)
formatDec(Ten4, displaydigits=4)
formatDec(Ten4, scientific=TRUE)
formatDec(Ten4, scientific=TRUE, displaydigits=4)


