\name{pochMpfr}
\alias{pochMpfr}
\title{Pochhammer Symbol aka Rising / Ascending Factorial}
\description{
  Compute the Pochhammer symbol or \dQuote{rising factorial}, also called
  the \dQuote{Pochhammer function}, \dQuote{Pochhammer polynomial},
  \dQuote{ascending factorial}, \dQuote{rising sequential product} or
  \dQuote{upper factorial},

  \deqn{x^{(n)}=x(x+1)(x+2)\cdots(x+n-1)= \frac{(x+n-1)!}{(x-1)!} = \frac{\Gamma(x+n)}{\Gamma(x)}.
  }{x^(n) = x(x+1)(x+2)...(x+n-1) = (x+n-1)! / (x-1)! = Gamma(x+n) / Gamma(x).}

}
\usage{
pochMpfr(a, n)
}
\arguments{
  \item{a}{a numeric or \code{\linkS4class{mpfr}} vector.}
  \item{n}{an integer (vector).}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  an \code{\linkS4class{mpfr}} vector as \code{a}.
}
\seealso{
  \code{\link{choose}(n,m)} computes the binomial coefficient
  \eqn{C_{n,m}}{C[n,m]} which can also be expressed via Pochhammer
  symbol as
  \eqn{C_{n,m} = (n-m+1)^{(m)}/m!}{C[n,m] = (n-m+1)^(m) / m!}.

  \code{\link{factorialMpfr}}.
}
\examples{
pochMpfr(100, 4) == 100*101*102*103 # TRUE
a <- 100:110
pochMpfr(a, 10) # exact (but too high precision)
x <- mpfr(a, 70)# should be enough
pochMpfr(x, 10) # the same as above (needing only 70 bits)
stopifnot(pochMpfr(a, 10) == pochMpfr(x, 10))
}
\keyword{arith}
