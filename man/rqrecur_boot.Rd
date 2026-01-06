\name{rqrecur_boot}
\alias{rqrecur_boot}
\title{Bootstrap a fitted recurrent-episode quantile regression model}
\description{
This function performs a \strong{cluster (subject-level)} bootstrap: subjects are resampled with
replacement and their episode rows are carried along.
}
\details{
The bootstrap is used to form \strong{pointwise} confidence intervals for the coefficient
functions \eqn{\beta(\tau)} across the fitted \eqn{\tau} grid.
}
\usage{
rqrecur_boot(object, data, B = 200, seed = NULL, progress = TRUE, max_tries = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{"rqrecur"} object.}
\item{data}{The original long-format data used for fitting.}
\item{B}{Number of bootstrap replicates.}
\item{seed}{Optional integer seed for reproducibility.}
\item{progress}{Logical; if \code{TRUE}, show a progress bar.}
\item{max_tries}{Maximum number of attempted refits (to recover from occasional bootstrap failures).
Default is \code{5 * B}.}
\item{...}{Unused.}
}
\value{
The input \code{object}, augmented with a \code{bootstrap} element.
}
\seealso{\code{\link{rqrecur}}, \code{\link[stats]{confint}}}
