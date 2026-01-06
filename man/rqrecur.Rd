\name{rqrecur}
\alias{rqrecur}
\title{Fit the quantile regression model for recurrent episode lengths}
\description{
Implements the two-step estimation procedure described in
\emph{Exploring the Heterogeneity in Recurrent Episode Lengths Based On Quantile Regression}
(Liu, Umpierrez, and Peng).
}
\details{
The method targets the conditional quantiles of episode lengths:
\deqn{Q_{X_{ij}}(\tau \mid Z_i, L_{ij}) = \beta_1(\tau)^\top Z_i + \beta_2(\tau)^\top L_{ij}}

and addresses dependent truncation, dependent censoring, and informative cluster size
via a two-step strategy:
\itemize{
  \item estimate \eqn{\beta_2(\tau)} from within-subject differences; and
  \item estimate \eqn{\beta_1(\tau)} by weighted quantile regression.
}
}
\usage{
rqrecur(
  formula,
  data,
  id,
  tau = c(0.25, 0.5, 0.75),
  td_vars,
  episode = NULL,
  complete = TRUE,
  status = NULL,
  status_complete = 1,
  censor_time = NULL,
  w1 = NULL,
  eps = 1e-8,
  boot = FALSE,
  B = 200,
  bootstrap = NULL,
  seed = NULL,
  progress = TRUE,
  ...
)
}
\arguments{
\item{formula}{A formula of the form \code{episode_length ~ x1 + x2 + ...} including both
time-independent and time-dependent covariates. Currently only additive terms are supported
(no interactions or transformations).}

\item{data}{A data.frame in long format with one row per (observed) episode.}

\item{id}{Name of the subject id column.}

\item{tau}{Quantile level(s) in \code{(0, 1)}. A numeric vector.}

\item{td_vars}{Character vector giving the \emph{time-dependent} covariate names among the predictors
in \code{formula}. The remaining predictors are treated as time-independent.
Use \code{character(0)} if there are no time-dependent covariates.}

\item{episode}{Optional name of an episode-order column. If \code{NULL}, the function will use row order within subject.}

\item{complete}{Logical. If \code{TRUE}, \code{data} is assumed to contain only complete episodes.
If \code{FALSE}, incomplete (censored) episodes are removed using \code{status}.}

\item{status}{Episode censoring indicator used when \code{complete = FALSE}. Either a column name
in \code{data} or a logical/numeric vector of length \code{nrow(data)}. For a numeric
indicator, \code{status_complete} specifies which value corresponds to a complete episode.}

\item{status_complete}{Value in \code{status} that indicates a complete episode when \code{status} is
numeric. Default is \code{1} (1=complete, 0=censored).}

\item{censor_time}{Censoring/follow-up times.
\itemize{
  \item If a numeric vector, it is treated as the cohort censoring times \eqn{C_1,\dots,C_n}.
  \item If a string, it is interpreted as a column in \code{data}; the function will extract it from \code{data}.
  \item If \code{NULL}, the function will look for a \code{censor_time} column in \code{data}.
}}

\item{w1}{Name of the column containing \eqn{W_{i1}} (time from study entry to start of
the first episode). If \code{NULL}, the function tries \code{"T1"} then \code{"w1"}, then uses the \code{"end"}
time of the first episode if available. If no truncation correction is performed, \code{w1} is not required.}

\item{eps}{Small positive constant used to avoid division by zero in \code{1/G(t)}. Default \code{1e-8}.}

\item{boot}{Logical. If \code{TRUE}, run a \strong{subject-level} bootstrap inside \code{rqrecur()} and attach
pointwise confidence intervals to the returned object.}

\item{B}{Number of bootstrap replicates (only used when \code{boot = TRUE}). Default is \code{200}.}

\item{bootstrap}{Deprecated alias for \code{B}. If provided, it overrides \code{B} and also sets
\code{boot = (bootstrap > 0)} for backward compatibility.}

\item{seed}{Optional integer seed for reproducible bootstrapping.}

\item{progress}{Logical; if \code{TRUE}, show a progress bar during bootstrapping.}

\item{...}{Currently unused (reserved).}
}
\value{
An object of class \code{"rqrecur"} containing coefficient estimates for each \code{tau}.
}
\examples{
dat <- read_example_sim_data()
fit <- rqrecur(
  epi_length ~ z1 + z2 + z3,
  data = dat,
  id = "id",
  td_vars = "z3",
  tau = c(0.25, 0.5, 0.75),
  # In real use, provide censor_time from your full cohort
  censor_time = NULL
)

# Predict quantiles for new covariates
newdat <- data.frame(z1 = 1.3, z2 = 1, z3 = 5)
predict(fit, newdata = newdat)
}
