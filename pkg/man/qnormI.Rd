\name{qnormI}
\alias{qnormI}
\title{Gaussian / Normal Quantiles \code{qnorm()} via Inversion}
\description{
  Compute Gaussian or Normal Quantiles \code{\link{qnorm}(p, *)} via
  inversion of our \dQuote{mpfr-ified} arbitrary accurate
  \code{\link[Rmpfr]{pnorm}()}, using our \code{\link{unirootR}()} root
  finder.
}
\usage{
qnormI(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE,
       trace = 0, verbose = as.logical(trace),
       tol,
       give.full = FALSE,
       \dots)
}
\arguments{% ~/R/D/r-devel/R/src/library/stats/man/Normal.Rd <<<<<
  \item{p}{vector of probabilities.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}.}
  %%--
  \item{trace}{integer passed to \code{\link{unirootR}()}.  If positive,
    information about a search interval extension will be printed to the console.}
  \item{verbose}{logical indicating if progress details should be printed
    to the console.}
  \item{tol}{the desired accuracy (convergence tolerance).}
  \item{give.full}{logical indicating if the \emph{full} result of
    \code{\link{unirootR}()} should be returned (when applicable).}
  \item{\dots}{optional further arguments passed to \code{\link{unirootR}()}
    such as \code{maxiter}, \code{verbDigits}, \code{check.conv},
    \code{warn.no.convergence}, and \code{epsC}.}
}
\value{
  If \code{give.full} is true, return a \code{\link{list}}, say \code{r}, of
  \code{\link{unirootR}(.)} results, with \code{length(r) == length(p)}.

  Otherwise, return a \dQuote{numeric vector} like \code{p}, e.g., of
  \code{class "mpfr"} when \code{p} is.
}
\author{Martin Maechler}

\seealso{
 Standard \R's \code{\link[stats]{qnorm}}.
}
\examples{
p  <- (0:32)/32
lp <- -c(1000, 500, 200, 100, 50, 20:1, 2^-(1:8))

for(logp  in c(FALSE,TRUE)) {
  pp <- if(logp) lp else p
  mp <- mpfr(pp, precBits = 80) # precBits = 128 gave "the same"
  for(l.tail in c(FALSE,TRUE)) {
      qn <- qnorm (pp, lower.tail = l.tail, log.p = logp)
     qnI <- qnormI(pp, lower.tail = l.tail, log.p = logp, tol = 2.3e-16)
     qnM <- qnormI(mp, lower.tail = l.tail, log.p = logp, tol =   1e-20)
     cat(sprintf("Accuracy of qnorm(*, lower.t=\%-5s, log.p=\%-5s): \%s\n",
                 l.tail, logp, format(all.equal(qnM, qn, tol=0))))
     stopifnot(exprs = {
        all.equal(qn,  qnI, tol = if(logp) 4e-14 else 8e-16)
        all.equal(qnM, qnI, tol = 2.3e-16)
     })
  }
}

## Boundary -- from limits in mpfr maximal exponent range!
## 1) Use maximal ranges:
old_eranges <- .mpfr_erange() # typically -/+ 2^30
(myERng <- (1-2^-52) * .mpfr_erange(c("min.emin","max.emax")))
doIncr <- !isTRUE(all.equal(myERng, old_eranges))
if(doIncr) .mpfr_erange_set(value = myERng)

## The largest quantile for which our mpfr-ized qnorm() does *NOT* underflow :
(qM <- mpfr("2528468770.343293436810768159197281514373932815851856314908753969469064"))
##           1 3 5 7 9  1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1 3 5 7 9 1
##                     10        20        30        40        50        60        70        80
(pM <- pnorm(-qM)) # precision  223   bits
 ## 7.64890682545699845135633468495894619457903458325606933043966616334460003e-1388255822130839040
log(pM) # ... 233 bits: -3196577161300663205.85759196211156141481203239336338270527868730785529034
%% FIXME : here, the internal qnInt() gives (-Inf, -Inf) :
try( qnormI(pM) ) ## Error: lower < upper not fulfilled (evt. TODO)
## but this works
(qnI <- qnormI(log(pM), log.p=TRUE))
all.equal(-qM, qnI, tol = 0) # << show how close
stopifnot( all.equal(-qM, qnI, tol = 1e-18) )# see 1.084202e-19

if(FALSE) # this (*SLOW*) gives 21 x the exact *same* result --- FIXME!
qnormI(log(pM) * (2:22), log.p=TRUE)
%% 21 'mpfr' numbers of precision  233   bits
%%  [1] -2528468770.34329343681076815919728151437393281585185631490875396946906388
%%  [2] -2528468770.34329343681076815919728151437393281585185631490875396946906388
%%  ...............
%% [21] -2528468770.34329343681076815919728151437393281585185631490875396946906388


## reset to previous status if needed
if(doIncr) .mpfr_erange_set( , old_eranges)
}
\keyword{distribution}
\keyword{math}
