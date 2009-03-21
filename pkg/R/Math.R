#### Define mpfr methods for Math  group functions
####                        ======

### "Arith", "Compare",..., are in ./Arith.R
###  ----                            ~~~~~~~

## [1] "abs"    "sign"    "sqrt"    "ceiling" "floor" "trunc" "cummax"
## [8] "cummin" "cumprod" "cumsum"  "exp"     "expm1" "log"   "log10"
##[15] "log2"   "log1p"   "cos"     "cosh"    "sin"   "sinh"  "tan"
##[22] "tanh"   "acos"    "acosh"   "asin"    "asinh" "atan"  "atanh"
##[29] "gamma"  "lgamma"  "digamma" "trigamma"

if(FALSE) ## here are the individual function
    dput(getGroupMembers("Math"))

## Uniform interface to C:
##
## Pass integer code to call and do the rest in C
## Codes from ~/R/D/r-devel/R/src/main/names.c :
.Math.codes <-
    c(
      "floor" =     1,
      "ceiling" =   2,
      "sqrt" =      3,
      "sign" =      4,
      "exp" =      10,
      "expm1" =    11,
      "log1p" =    12,
      "cos" =      20,
      "sin" =      21,
      "tan" =      22,
      "acos" =     23,
      "asin" =     24,
      "cosh" =     30,
      "sinh" =     31,
      "tanh" =     32,
      "acosh" =    33,
      "asinh" =    34,
      "atanh" =    35,
      "lgamma" =   40,
      "gamma" =    41,
      "digamma" =  42,
      "trigamma" = 43)

.Math.gen <- getGroupMembers("Math")

## Those "Math" group generics that are not in the do_math1 table above
if(FALSE)
.Math.gen[!(.Math.gen %in% names(.Math.codes))] # grouped "semantically":
## "abs" "atan"   {no problem}
## "trunc"        has ... for methods
## "log"
## "log10" "log2" {no problem}
##
## "cummax" "cummin" "cumprod" "cumsum"  { no problem; but  related to
##                                         max,min,..: --> "Summary" group }

.Math.codes <-
    c(.Math.codes,
      "trunc" = 0, "atan" = 25, # "abs" has own method!
      "log" = 13, "log2" = 14, "log10" = 15,
      "cummax" = 51, "cummin" = 52, "cumprod" = 53, "cumsum" = 54,
      ## These are *NOT* in R's  Math group, but 1-argument math functions
      ## available in the mpfr - library:
      "erf" = 101, "erfc" = 102, "zeta" = 104, "Eint" = 106,
      "j0" = 111, "j1" = 112, "y0" = 113, "y1" = 114)
storage.mode(.Math.codes) <- "integer"


## A few ones have a very simple method:
setMethod("sign", "mpfr",
	  function(x) sapply(x, function(e) e@sign))

setMethod("abs", "mpfr",
	  function(x) {
	      for(i in seq_along(x)) x[[i]]@sign <- 1L
	      x
	  })

## "log" is still special with its 'base' :
setMethod("log", signature(x = "mpfr"),
	  function(x, base) {
	      if(!missing(base) && base != exp(1))
		  stop("base != exp(1) is not yet implemented")
	      new("mpfr", .Call("Math_mpfr", x, .Math.codes["log"],
				PACKAGE="Rmpfr"))
	  })

setMethod("Math", signature(x = "mpfr"),
	  function(x) {
	      new("mpfr", .Call("Math_mpfr", x, .Math.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })
