\name{mpfr}
\alias{mpfr}
\alias{Const}
\title{Create "mpfr" Numbers (Objects)}
\description{
  Create multiple (i.e. typically \emph{high}) precision numbers, to be
  used in arithmetic and mathematical computations with \R.
}
\usage{
mpfr(x, precBits)
Const(name = c("pi", "gamma", "catalan"), prec = 120L)
}
\arguments{
  \item{x}{a \code{\link{numeric}} vector.}
  \item{precBits, prec}{a number, the maximal precision to be used, in
    \bold{\emph{bits}}; i.e. \code{53} corresponds to double precision.
    Must be at least 2.}
  \item{name}{a string specifying the mpfrlib - internal constant
    computation.  \code{"gamma"} is Euler's gamma \eqn{(\gamma)}{}, and
    \code{"catalan"} Catalan's constant.}
}
\value{
  an object of (S4) class \code{\linkS4class{mpfr}} which the user
  should just as a normal numeric vector.
}
\author{Martin Maechler}
\seealso{The class documentation \code{\linkS4class{mpfr}} contains more
  details.
}
\examples{
mpfr(pi, 120) ## the double-precision pi "translated" to 120-bit precision

pi. <- Const("pi", prec = 260) # pi "computed" to correct 260-bit precision
pi. # nicely prints 80 digits [260 * log10(2) ~= 78.3 ~ 80]

Const("gamma",   128L) # 0.5772...
Const("catalan", 128L) # 0.9159...

x <- mpfr(0:7, 100)/7 # a more precise version of  k/7, k=0,..,7
x
1 / x
}
\keyword{classes}

