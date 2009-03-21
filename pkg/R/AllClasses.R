#### All Class Definitions in package  "Rmpfr"

### NB:	 Use   /usr/local/app/R/R_local/src/Brobdingnag/R/brob.R
###					    -----------
### as a partial role image

setClass("mpfr1", ## a single Multi-precision float number
	 representation(prec = "integer", # precision in bits
			exp = "integer",  # exponent
			sign= "integer",  # signum
			d = "integer"),	  # the mantissa as a vector of long
	 validity = function(object) {
	     if(length(pr <- object@prec) != 1 || is.na(pr) || pr < 2)
		 "invalid 'prec' slot"
	     else if(length(ex <- object@prec) != 1)
		 "invalid 'exp' slot"
	     else if(length(sig <- object@sign) != 1 || is.na(sig) ||
		     abs(sig) > 1)
		 "invalid 'sign' slot"
	     else if(length(d <- object@d) != ceiling(pr / 32))
		 "length('d' slot) does not match 'prec'"
	     else TRUE
	 })

setClass("mpfr", ## a *vector* of "mpfr1", i.e., multi-precision float numbers
	 contains = "list", ## of "mpfr1" entries:
	 validity = function(object) {
	     if(all(sapply(object, class) == "mpfr1"))
		 return(TRUE)
	     ## else
		 "Not all components are of class 'mpfr1'"
	 })
