\name{factorialMpfr}
\alias{factorialMpfr}
\title{Factorial 'n!'  in Arbitrary Precision}
\description{
  Efficiently compute \eqn{n!} in arbitrary precision,
  using the MPFR-internal implementation.
  This is mathematically (but not numerically) the same as
  \eqn{\Gamma(n+1)}{Gamma(n+1)}.
}
\usage{
factorialMpfr(n, precBits = ceiling(lgamma(n+1)/log(2)))
}
\arguments{
  \item{n}{non-negative integer.}
  \item{precBits}{desired precision in bits (\dQuote{binary digits});
    the default sets the precision high enough for the result to be \emph{exact}.}
}
\value{
  a number of (S4) class \code{\linkS4class{mpfr}}.
}
\seealso{
  \code{\link{factorial}} and \code{\link{gamma}} in base \R;
  \code{\link{pochMpfr}}.
}
\examples{
factorialMpfr(200)

f1000 <- factorialMpfr(1000)
stopifnot(1e-15 > abs(as.numeric(1 - lfactorial(1000)/log(f1000))))

## >> FIXME: we have a factorial(mpfr(...)) example at the end of
## >> of ./mpfr-class.Rd  ... that is now "sub-optimal"
## HOWEVER: Measurements show only *small* efficiency gain of ~ 10\% :
system.time(replicate(8, f1e4 <- factorialMpfr(10000)))
system.time(replicate(8, f.1e4 <- factorial(mpfr(10000,
                            prec=1+lfactorial(10000)/log(2)))))
}
\keyword{arith}
