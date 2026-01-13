\name{predict.rqrecur}
\alias{predict.rqrecur}
\title{Predict episode-length quantiles for new covariates}
\usage{
\method{predict}{rqrecur}(object, newdata, tau = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{"rqrecur"} object.}
\item{newdata}{A data.frame containing the covariates used in the fit.}
\item{tau}{Optional subset of fitted quantiles to predict.}
\item{...}{Unused.}
}
\value{
A numeric matrix with one row per row of \code{newdata} and one column per \code{tau}.
}
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

newdat <- data.frame(z1 = 1.3, z2 = 1, z3 = 5)
predict(fit, newdata = newdat)
}
