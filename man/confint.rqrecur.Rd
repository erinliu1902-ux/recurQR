\name{confint.rqrecur}
\alias{confint.rqrecur}
\title{Pointwise confidence intervals for coefficient functions}
\description{
Computes pointwise bootstrap confidence intervals for \eqn{\beta(\tau)}.
}
\arguments{
\item{object}{A fitted \code{"rqrecur"} object with bootstrap results (see \code{\link{rqrecur_boot}} or
fit with \code{boot = TRUE} in \code{\link{rqrecur}}).}
\item{parm}{Optional coefficient name(s). If \code{NULL}, returns all coefficients.}
\item{level}{Confidence level (default 0.95).}
\item{...}{Unused.}
}
\usage{
\method{confint}{rqrecur}(object, parm = NULL, level = 0.95, ...)
}
\value{
A list with elements \code{estimate}, \code{lower}, \code{upper}, and \code{level}.
}
\seealso{\code{\link{rqrecur}}, \code{\link{rqrecur_boot}}}
