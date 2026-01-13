\docType{package}
\name{recurQR-package}
\alias{recurQR-package}
\alias{recurQR}
\title{recurQR: Quantile Regression for Recurrent Episode Lengths}
\description{
The \pkg{recurQR} package implements the two-step estimation procedure proposed in
Liu, Umpierrez, and Peng (2025) for modeling recurrent episode lengths via quantile regression
while addressing dependent truncation, dependent censoring, and informative cluster size.
}
\details{
The main modeling function is \code{\link{rqrecur}}.

Other helpful functions/methods include:
\itemize{
  \item \code{\link{read_example_sim_data}} (load built-in example data),
  \item \code{\link{rqrecur_boot}} (subject-level bootstrap for pointwise CIs),
  \item \code{\link[stats]{predict}} and \code{\link{predict.rqrecur}},
  \item \code{\link[graphics]{plot}} and \code{\link{plot.rqrecur}},
  \item \code{\link[stats]{confint}} and \code{\link{confint.rqrecur}}.
}

A worked example is available in the package vignette:
\code{vignette("recurQR")}.
}
\author{
Yi Liu; Guillermo E. Umpierrez; Limin Peng \email{lpeng@emory.edu}
}
\keyword{package}
\examples{
dat <- read_example_sim_data()
fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = 0.5,
  censor_time = NULL
)
coef(fit)
}
