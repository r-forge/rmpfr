#### Define mpfr methods for Arith + Compare + Logic	group functions
####			    ======   =======   =====

### "Math" are done in ./Math.R ,  "Summary" in ./Summary.R
###  ----		 ~~~~~~     -------       ~~~~~~~~~
### NB:	 Look at /usr/local/app/R/R_local/src/Brobdingnag/R/brob.R
###					      -----------

if(FALSE) ## here are the individual function
getGroupMembers("Ops")
if(FALSE) ## here are the individual function
str(list(Arith	 = getGroupMembers("Arith"),
	 Compare = getGroupMembers("Compare"),
	 Logic	 = getGroupMembers("Logic")), vec.len = 20)

setMethod("Ops", signature(e1 = "mpfr", e2 = "ANY"),
	  function(e1, e2) callGeneric(e1, as(e2, "numeric")))
setMethod("Ops", signature(e1 = "ANY", e2 = "mpfr"),
	  function(e1, e2) callGeneric(as(e1, "numeric"), e2))

setMethod("Logic", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) callGeneric(as(e1, "numeric"),
				       as(e2, "numeric")))

###-- 2) ----------- Arith --------------------------------------------------

.mpfr.negative <- function(x) {
    for(i in seq_along(x)) x[[i]]@sign <- -x[[i]]@sign
    x
}

setMethod("Arith",signature(e1 = "mpfr", e2="missing"),
	  function(e1,e2) {
	    switch(.Generic,
		   "+" = e1,
		   "-" = .mpfr.negative(e1),
		   stop(paste("Unary operator", .Generic,
			      "not defined for \"mpfr\" numbers"))
		   )
	  } )


.Arith.codes <-
    c("+" = 1, "-" = 2, "*" = 3, "^" = 4, "%%" = 5, "%/%" =6, "/" = 7)
storage.mode(.Arith.codes) <- "integer"

setMethod("Arith", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) {
	      new("mpfr", .Call("Arith_mpfr", e1, e2, .Arith.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })

setMethod("Arith", signature(e1 = "mpfr", e2 = "integer"),
	  function(e1, e2) {
	      new("mpfr", .Call("Arith_mpfr_i", e1, e2, .Arith.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })
setMethod("Arith", signature(e1 = "integer", e2 = "mpfr"),
	  function(e1, e2) {
	      new("mpfr", .Call("Arith_i_mpfr", e1, e2, .Arith.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })

setMethod("Arith", signature(e1 = "mpfr", e2 = "numeric"),# not "integer"
	  function(e1, e2) {
	      new("mpfr", .Call("Arith_mpfr_d", e1, e2, .Arith.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })
setMethod("Arith", signature(e1 = "numeric", e2 = "mpfr"),# not "integer
	  function(e1, e2) {
	      new("mpfr", .Call("Arith_d_mpfr", e1, e2, .Arith.codes[.Generic],
				PACKAGE="Rmpfr"))
	  })


###-- 3) ----------- Compare --------------------------------------------------

.Compare.codes <- c("==" = 1, ">" = 2, "<" = 3, "!=" = 4, "<=" = 5, ">=" =6)
storage.mode(.Compare.codes) <- "integer"

setMethod("Compare", signature(e1 = "mpfr", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call("Compare_mpfr", e1, e2, .Compare.codes[.Generic],
		    PACKAGE="Rmpfr")
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 = "integer"),
	  function(e1, e2) {
	      .Call("Compare_mpfr_i", e1, e2, .Compare.codes[.Generic],
		    PACKAGE="Rmpfr")
	  })

setMethod("Compare", signature(e1 = "mpfr", e2 =  "numeric"),# not "integer"
	  function(e1, e2) {
	      .Call("Compare_mpfr_d", e1, e2, .Compare.codes[.Generic],
		    PACKAGE="Rmpfr")
	  })

setMethod("Compare", signature(e1 = "integer", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call("Compare_mpfr_i", e2, e1,
		    .Compare.codes[c(1, 3:2, 4, 6:5)][.Generic],
		    PACKAGE="Rmpfr")
	  })

setMethod("Compare", signature(e1 = "numeric", e2 = "mpfr"),
	  function(e1, e2) {
	      .Call("Compare_mpfr_d", e2, e1,
		    .Compare.codes[c(1, 3:2, 4, 6:5)][.Generic],
		    PACKAGE="Rmpfr")
	  })



## ".Brob.arith" <- function(e1,e2){
##   switch(.Generic,
##          "+" = .Brob.add  (e1, e2),
##          "-" = .Brob.add  (e1, .Brob.negative(as.brob(e2))),
##          "*" = .Brob.mult (e1, e2),
##          "/" = .Brob.mult (e1, .Brob.inverse(as.brob(e2))),
##          "^" = .Brob.power(e1, e2),
##          stop(paste("binary operator \"", .Generic, "\" not defined for Brobdingnagian numbers"))
##          ) }

## setMethod("Arith", signature(e1 = "brob", e2="ANY"), .Brob.arith)
## setMethod("Arith", signature(e1 = "ANY", e2="brob"), .Brob.arith)
## setMethod("Arith", signature(e1 = "brob", e2="brob"), .Brob.arith)
