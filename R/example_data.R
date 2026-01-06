#' Read the example simulation dataset shipped with the package
#'
#' The package includes a small example dataset (an `.rds` file) derived from the authors'
#' simulation code. The dataset is in long format with one row per episode.
#'
#' @return A data.frame.
#' @export
read_example_sim_data <- function() {
  path <- system.file("extdata", "example_sim_data.rds", package = "recurQR")
  if (path == "") stop("example_sim_data.rds not found inside the installed package.", call. = FALSE)
  readRDS(path)
}
