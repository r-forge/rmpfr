% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Rmpfr/inst/NEWS.Rd"))
\name{NEWS}
\title{Rmpfr News}
\encoding{UTF-8}

\section{Changes in version 0.6-0 [2015-07-04]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{str()} is now based on \code{formatMpfr()} and nicely
      shows numbers also out of the double precision range.
      \item \code{matrix(mp, ..)} now also works when \code{mp} is of
      class \code{"mpfr"}.
      \item new matrix \code{norm()} for several \code{kind}s.
      \item new functions \code{formatHex()} and \code{formatBin()}
      thanks to Rich Heiberger.
      \item \code{mpfr(x)} for character \code{x} now auto-detects
      binary or hexadecimal numbers and hence acts as \emph{inverse} of
      \code{formatBin} and \code{formatHex}.
      \item \code{str(<mpfr>, ....)} now chooses a smart default for
      optional argument \code{vec.len}.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{pnorm(.)}, \code{j0()} and similar special functions
      now preserve \code{mpfrMatrix} and \code{mpfrArray} classes.

      \item similarly, \code{is.finite()} etc keep the
      \code{dim()}ensionality for \code{"mpfrArray"} arguments.

      \item \code{mpfr("0xabc", base=16)} and \code{mpfr("0b101",
	base=2)} and corresponding \code{getPrec()} now give the correct
      precBits instead of (8 or 2, respectively) too many.

      \item \code{str(<0-length mpfr>)} now works correctly.
    }
  }
}



\section{Changes in version 0.5-7 [2014-11-27]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item \code{as.integer()} now rounds \dQuote{to zero} as for
      regular \R numbers (it accidentally did round \dQuote{to nearest}
      previously).
    }
  }
}

\section{Changes in version 0.5-6 [2014-09-05]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item experimental \code{mpfrImport()}, \code{mpfrXport()}
      utilities -- as we found cases, where save() \code{"mpfr"} objects
      were \emph{not} portable between different platforms.

      \item \code{as(*,"mpfr")} now also supports rounding mode
      \code{"A"} (\dQuote{\bold{A}way from zero}).

      \item Several hidden low level utilities also get a
      \code{rnd.mode} option.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}

\section{Changes in version 0.5-5 [2014-06-19, r190]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item The result of \code{integrateR()} now prints even if a
      warning happened.

      \item \code{pbetaI(x, a,b)}, the arbitrarily accurate
      \code{pbeta()} computation for \emph{integer} \eqn{a} and \eqn{b},
      now works for larger \eqn{(a,b)}.

      \item Newly providing \code{mpfr}-versions of \code{dbinom()},
      \code{dpois()}, and \code{dnorm()}.

      \item New utility functions \code{mpfr_default_prec()},
      \code{.mpfr.minPrec()}, etc, to get, check, set default exponent
      ranges and precision.

      \item New \code{sinpi()}, \code{cospi()} etc, notably for \R >= 3.0.1.
    }
  }
  %% \subsection{BUG FIXES}{
  %%   \itemize{
  %%     \item .
  %%   }
  %% }
}

\section{Changes in version 0.5-4 [2013-10-22, r173]}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item .
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item .
    }
  }
}
