\name{mpfrMatrix}
\title{Classes "mpfrMatrix" and "mpfrArray"}
\docType{class}
\alias{mpfrMatrix-class}
\alias{mpfrArray-class}
%
\description{The classes \code{"mpfrMatrix"} and \code{"mpfrArray"} are,
  analogously to the \pkg{base} \code{\link{matrix}} and \code{\link{array}}
  functions and classes simply \dQuote{numbers} of class
  \code{\linkS4class{mpfr}} with an additional \code{Dim} and
  \code{Dimnames} slot.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("mpfrMatrix", ...)} or \code{new("mpfrArray", ...)}, or also
  by \code{dim(x) <- dd} where \code{x} is a an
  \code{\linkS4class{mpfr}} \dQuote{number vector}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{as for the \code{\linkS4class{mpfr}} class, a
      \code{"list"} of  \code{\linkS4class{mpfr1}} numbers.}
    \item{\code{Dim}:}{of class \code{"integer"}, specifying the array
      dimension.}
    \item{\code{DimNames}:}{of class \code{"list"} and the same length
      as \code{Dim}, each list component either \code{\link{NULL}} or a
      \code{\link{character}} vector of length \code{Dim[j]}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{mpfrArray}"}, directly.
Class \code{"\linkS4class{mpfr}"}, by class "mpfrArray", distance 2.
Class \code{"\linkS4class{list}"}, by class "mpfrArray", distance 3.
Class \code{"\linkS4class{vector}"}, by class "mpfrArray", distance 4.
}
%% FIXME in time
% \section{Methods}{
% No methods defined with class "mpfrMatrix" in the signature.
% }
\author{Martin Maechler}
% \note{ ~~further notes~~ }
% \seealso{
% 	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
% 	or \code{\linkS4class{CLASSNAME}} for links to other classes
% }
\examples{
showClass("mpfrMatrix")

validObject(mm <- new("mpfrMatrix"))
validObject(aa <- new("mpfrArray"))

m6 <- new("mpfrMatrix", mpfr(1:6, 128), Dim = c(2L, 3L))
validObject(m6)
m6
which(m6 == 3, arr.ind = TRUE) # |-->  (1, 2)
}
\keyword{classes}
