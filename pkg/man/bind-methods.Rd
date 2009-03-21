\name{cbind-methods}
\title{"mpfr" Methods for Functions cbind() \& rbind()}}
\docType{methods}
\alias{cbind-methods}
\alias{rbind-methods}
\alias{cbind,ANY-method}
\alias{cbind,Mnumber-method}
\alias{rbind,ANY-method}
\alias{rbind,Mnumber-method}
\description{
  \code{\link{cbind}} and \code{\link{rbind}} methods for signature
  \code{...} (see \code{\link[methods]{dotsMethods}}
  are provided for class \code{\linkS4class{Mnumber}}, i.e., for binding
  numeric vectors and class \code{"\linkS4class{mpfr}"} vectors and
  matrices (\code{"\linkS4class{mpfrMatrix}"}) together.
}
\section{Methods}{
  \describe{
    \item{... = "Mnumber"}{is used to (c|r)bind multiprecision \dQuote{numbers}
      (inheriting from class \code{"\linkS4class{mpfr}"})
      together, maybe combined with simple numeric vectors.}
    \item{... = "ANY"}{reverts to \code{\link[base]{cbind}} and \code{rbind}
      from package \pkg{base}.}
  }
}
\keyword{methods}

