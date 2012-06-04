\name{optimizeR}
\alias{optimizeR}
\title{High Precisione One-Dimensional Optimization}
%% Compare with ~/R/D/r-devel/R/src/library/stats/man/optimize.Rd
\description{
  \code{optimizeR} searches the intervalfrom
  \code{lower} to \code{upper} for a minimum % TODO ? or maximum
  of the function \code{f} with respect to its first argument.
}
\usage{
optimizeR(f, lower, upper, ..., tol = 1e-20,
          method = c("GoldenRatio"),
          extraBits = 30, precBits = -log2(tol) + extraBits, maxiter = 1000)
}
\arguments{
  \item{f}{the function to be minimized. \code{f(x)} must work
    \dQuote{in \pkg{Rmpfr} arithmetic} for \code{optimizer()} to make sense.}
    % optimized .. The function is
    % either minimized or maximized over its first argument
    % depending on the value of \code{maximum}.}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
  \item{tol}{the desired accuracy, typically higher than double
    precision, i.e., \code{tol < 2e-16}.}
  \item{method}{\code{\link{character}} string specifying the
    optimization method.  Currently, only "GoldenRatio" is available.}
  \item{extraBits}{only for default \code{precBits} construction: number
    of extra bits in addition those directly needed for \code{tol}.}
  \item{precBits}{number of bits to be used for
    \code{\link[Rmpfr]{mpfr}} numbers used internally.}
  \item{maxiter}{maximal number of iterations to be used.}
}
\details{
  %%\describe{\item{\code{"GoldenRatio"}:}{..........}
  Currently, there is method \code{"GoldenRatio"}, which from a given
  interval containing the solution, constructs the next point in the
  golden ratio between the interval boundaries.
  %%}% end{describe}
}
\value{
  A \code{\link{list}} with components \code{minimum} % (or \code{maximum})
  and \code{objective} which give the location of the minimum % (or maximum)
  and the value of the function at that point;
  \code{iter} specifiying the number of iterations, the logical
  \code{convergence} indicating if the iterations converged and
  \code{estim.prec} which is an estimate of the final precision (in \eqn{x}).
}

% \references{
% }
\author{
  Based on Hans W Borchert's \code{\link[pracma]{golden_ratio}};
  modifications by Martin Maechler.
}
% \note{
% }
\seealso{
  \R's standard \code{\link{optimize}}; \pkg{Rmpfr}'s \code{\link{unirootR}}.
}
\examples{
## Testing on the sine curve with 40 correct digits:
sol <- optimizeR(sin, 2, 6, tol = 10^-40)
str(sol)

print(rbind(2*sol$min/3, Const("pi", 256)),
      digits=40)
}
\keyword{optimize}
