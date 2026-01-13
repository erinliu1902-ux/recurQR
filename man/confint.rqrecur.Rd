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
\examples{
\dontrun{
dat <- read_example_sim_data()
fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = seq(0.2, 0.8, by = 0.1),
  censor_time = NULL
)

# Attach bootstrap results (pointwise CIs)
fit_b <- rqrecur_boot(fit, data = dat, B = 200, seed = 1)

# Pointwise confidence intervals for all coefficients
confint(fit_b)

# Or subset specific coefficients
confint(fit_b, parm = c("z1", "z3"))
}
}
