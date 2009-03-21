\name{mpfr-class}
\docType{class}
\title{Class "mpfr" of Multiple Precision Floating Point Numbers}
%
\alias{mpfr-class}
\alias{mpfr1-class}
%
\alias{[,mpfr,ANY,missing,missing-method}
\alias{[<-,mpfr,ANY,missing,ANY-method}
\alias{[<-,mpfr,ANY,missing,mpfr-method}

\alias{as.numeric,mpfr-method}
\alias{as.integer,mpfr-method}
\alias{coerce,mpfr,character-method}
\alias{coerce,mpfr,numeric-method}
\alias{coerce,mpfr1,numeric-method}
\alias{coerce,mpfr1,mpfr-method}
\alias{coerce,integer,mpfr-method}
\alias{coerce,numeric,mpfr-method}
\alias{coerce,numeric,mpfr1-method}

\alias{Ops,mpfr,ANY-method}
\alias{Ops,ANY,mpfr-method}
\alias{Arith,mpfr,missing-method}
\alias{Arith,mpfr,mpfr-method}
\alias{Arith,mpfr,integer-method}
\alias{Arith,mpfr,numeric-method}
\alias{Arith,integer,mpfr-method}
\alias{Arith,numeric,mpfr-method}
\alias{Compare,mpfr,mpfr-method}
\alias{Compare,mpfr,integer-method}
\alias{Compare,mpfr,numeric-method}
\alias{Compare,integer,mpfr-method}
\alias{Compare,numeric,mpfr-method}
\alias{Logic,mpfr,mpfr-method}
\alias{Summary,mpfr-method}
\alias{Math,mpfr-method}
\alias{abs,mpfr-method}
\alias{log,mpfr-method}
\alias{sign,mpfr-method}

\alias{format,mpfr-method}
\alias{is.finite,mpfr-method}
\alias{is.infinite,mpfr-method}
\alias{is.na,mpfr-method}
\alias{is.nan,mpfr-method}
\alias{unique,mpfr,missing-method}
%
\alias{all.equal,mpfr,ANY-method}
\alias{all.equal,ANY,mpfr-method}
%
\alias{show,mpfr-method}
\alias{show,mpfr1-method}
\alias{print.mpfr}
\alias{print.mpfr1}
%
\description{
  \code{"mpfr"} is the class of \bold{M}ultiple \bold{P}recision
  \bold{F}loatingpoint numbers with \bold{R}eliable arithmetic.

  For the high-level user, \code{"mpfr"} objects should behave
  as standard \R's \code{\link{numeric}} \emph{vectors}, just with
  prespecified (typically high) precision.
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mpfr}(<number>, precBits)}.
}
\section{Slots}{
  Internally, \code{"mpfr"} objects just contain standard \R
  \code{\link{list}}s where each list element is of class
  \code{"mpfr1"}, representing \emph{one} MPFR number, in a structure
  with four slots, very much parallelizing the C \code{struc} in the
  \code{mpfr} C library to which the \pkg{Rmpfr} package interfaces.

  An object of class \code{"mpfr1"} has slots
  \describe{
    \item{\code{prec}:}{\code{"integer"} specifying the maxmimal
      precision in \bold{bits}.}
    \item{\code{exp}:}{\code{"integer"} specifying the base-\bold{2}
      exponent of the number.}
    \item{\code{sign}:}{\code{"integer"}, typically \code{-1} or
      \code{1}, specifying the sign (i.e. \code{\link{sign}(.)}) of the
      number.}
    \item{\code{d}:}{an \code{"integer"} vector (of 32-bit
      \dQuote{limbs}) which corresponds to the full mantissa of the
      number.}
  }% desc
}
% \section{Extends}{
% Class \code{"\linkS4class{list}"}, from data part.
% Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
% }
\section{Methods}{
  \describe{
    \item{abs}{\code{signature(x = "mpfr")}: ... }

    \item{Ops}{\code{signature(e1 = "mpfr", e2 = "ANY")}: ... }
    \item{Ops}{\code{signature(e1 = "ANY", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "missing")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Arith}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Arith}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Arith}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "integer")}: ... }
    \item{Compare}{\code{signature(e1 = "mpfr", e2 = "numeric")}: ... }
    \item{Compare}{\code{signature(e1 = "integer", e2 = "mpfr")}: ... }
    \item{Compare}{\code{signature(e1 = "numeric", e2 = "mpfr")}: ... }

    \item{Logic}{\code{signature(e1 = "mpfr", e2 = "mpfr")}: ... }

    \item{Summary}{\code{signature(x = "mpfr")}: ... }

    \item{Math}{\code{signature(x = "mpfr")}: ... }

    \item{as.numeric}{\code{signature(x = "mpfr")}: ... }
    \item{[}{\code{signature(x = "mpfr", i = "ANY", j = "missing", drop = "missing")}: ... }
    \item{format}{\code{signature(x = "mpfr")}: ... }
    \item{is.finite}{\code{signature(x = "mpfr")}: ... }
    \item{is.infinite}{\code{signature(x = "mpfr")}: ... }
    \item{is.na}{\code{signature(x = "mpfr")}: ... }
    \item{is.nan}{\code{signature(x = "mpfr")}: ... }
    \item{log}{\code{signature(x = "mpfr")}: ... }
    \item{show}{\code{signature(object = "mpfr")}: ... }
    \item{sign}{\code{signature(x = "mpfr")}: ... }

    \item{all.equal}{\code{signature(target = "mpfr", current = "ANY")}, and}
    \item{all.equal}{\code{signature(target = "ANY", current = "mpfr")}:
    methods for numerical (approximate) equality,
    \code{\link[base]{all.equal}} of multiple precision numbers.  Note
    that it may make sense for these to use \code{tolerance} (argument)
    values of considerably less than double precision machine precision.}

    \item{coerce}{\code{signature(from = "numeric", to = "mpfr")}: ... }
    \item{coerce}{\code{signature(from = "mpfr", to = "numeric")}: ... }
    \item{coerce}{\code{signature(from = "mpfr", to = "character")}: ... }
  }
}
% \references{ ~put references to the literature/web site here ~ }
\author{Martin Maechler}
% \note{ ~~further notes~~ }

%  ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% 	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
% 	or \code{\linkS4class{CLASSNAME}} for links to other classes
% }
\examples{
## 30 digit precision
str(x <- mpfr(c(2:3, pi), prec = 30 * log2(10)))
x^2
x[1] / x[2] # 0.66666... ~ 30 digits
}
\keyword{classes}
