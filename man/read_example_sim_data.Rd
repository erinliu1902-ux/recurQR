\name{read_example_sim_data}
\alias{read_example_sim_data}
\title{Read the example simulation dataset shipped with the package}
\description{
Reads \code{inst/extdata/example_sim_data.rds} from the installed package.
}
\usage{
read_example_sim_data()
}
\value{
A data.frame in long format with one row per episode.
}
\examples{
dat <- read_example_sim_data()
dim(dat)
head(dat)
}
