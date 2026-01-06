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
