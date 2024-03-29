\name{mpfr-special-functions}
\alias{zeta}
\alias{Ei}
\alias{Li2}
\alias{erf}
\alias{erfc}
\title{Special Mathematical Functions (MPFR)}
\description{
  Special Mathematical Functions, supported by the MPFR Library.

  Note that additionally, all the \code{\link{Math}} and
  \code{\link{Math2}} group member functions are \dQuote{mpfr-ified}, too;
  ditto, for many more standard \R functions.  See see the methods listed
  in \code{\linkS4class{mpfr}} (aka \code{?`\link{mpfr-class}`}). % ./mpfr-class.Rd
}
\usage{
zeta(x)
Ei(x)
Li2(x)

erf(x)
erfc(x)
}
\arguments{
  \item{x}{a \code{\link{numeric}} or \code{\linkS4class{mpfr}} vector.}
}
\details{
  \code{zeta(x)} computes Riemann's Zeta function
  \eqn{\zeta(x)}{zeta(x)} important in analytical number theory and
  related fields. The traditional definition is
  \deqn{\zeta(x) = \sum_{n=1}^\infty \frac{1}{n^x}.}{Zeta(x) = sum[n=1..Inf; 1/(n^x)].}

  \code{Ei(x)} computes the \bold{e}xponential integral,
  \deqn{\int_{-\infty}^{x} \frac{e^t}{t} \; dt.}{Integral(-Inf,x; e^t/t dt).}

  \code{Li2(x)} computes the dilogarithm,
  \deqn{\int_{0}^{x} \frac{-log(1-t)}{t} \; dt.}{Integral(0,x; -log(1-t)/t dt).}

  \code{erf(x)} and \code{erfc(x)} are the error, respectively
  \bold{c}omplementary error function which are both reparametrizations
  of \code{\link{pnorm}}, \code{erf(x) = 2*pnorm(sqrt(2)*x)} and
  \code{erfc(x) = 2* pnorm(sqrt(2)*x, lower=FALSE)},
  and hence \pkg{Rmpfr} provides its own version of \code{\link{pnorm}}.
}
\value{
  A vector of the same length as \code{x}, of class \code{\linkS4class{mpfr}}.
}
\seealso{\code{\link[stats:Normal]{pnorm}} in standard package \pkg{stats};
  the class description \code{\linkS4class{mpfr}} mentioning the
  generic arithmetic and mathematical functions (\code{sin}, \code{log},
  \dots, etc) for which \code{"mpfr"} methods are available.

  Note the (integer order, non modified) Bessel functions \eqn{j_0()},
  \eqn{y_n()}, etc, named \code{\link{j0}, \link{yn}} etc, and Airy
  function \eqn{Ai()} in \link{Bessel_mpfr}.
}
\examples{
curve(Ei,  0, 5, n=2001)

## As we now require (mpfrVersion() >= "2.4.0"):
curve(Li2,  0,    5, n=2001)
curve(Li2, -2,   13, n=2000); abline(h=0,v=0, lty=3)
curve(Li2, -200,400, n=2000); abline(h=0,v=0, lty=3)

curve(erf, -3,3, col = "red", ylim = c(-1,2))
curve(erfc, add = TRUE, col = "blue")
abline(h=0, v=0, lty=3)
legend(-3,1, c("erf(x)", "erfc(x)"), col = c("red","blue"), lty=1)
}
\keyword{math}
