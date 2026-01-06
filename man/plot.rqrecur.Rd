\name{plot.rqrecur}
\alias{plot.rqrecur}
\title{Plot fitted coefficient functions with optional pointwise confidence intervals}
\description{
This is a convenience method for plotting the estimated coefficient functions
\eqn{\beta(\tau)} across the fitted \eqn{\tau} grid.
}
\details{
Confidence intervals are obtained from subject-level bootstrap samples (see
\code{\link{rqrecur_boot}} or fit with \code{boot = TRUE} in \code{\link{rqrecur}}).
}
\usage{
\method{plot}{rqrecur}(x, parm = NULL, level = 0.95, ci = TRUE, ncol = NULL, ask = NULL, ...)
}
\arguments{
\item{x}{A fitted \code{rqrecur} object.}
\item{parm}{Coefficient name(s) to plot. If \code{NULL}, plots all coefficients.}
\item{level}{Confidence level for pointwise intervals (default 0.95).}
\item{ci}{Logical; whether to draw confidence intervals when available.}
\item{ncol}{Number of columns in the multi-panel layout when plotting multiple coefficients.
If \code{NULL}, a square-ish layout is chosen automatically.}
\item{ask}{Logical; if \code{TRUE}, pause between pages when many panels are drawn.}
\item{...}{Passed to the underlying \code{plot()} call (e.g., \code{ylim}, \code{lwd}).}
}
\value{
Coefficients plot.
}
\seealso{\code{\link{rqrecur}}, \code{\link{rqrecur_boot}}, \code{\link[stats]{confint}}}
