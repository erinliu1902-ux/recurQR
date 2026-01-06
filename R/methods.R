#' @export
print.rqrecur <- function(x, ...) {
  cat("<rqrecur> Two-step quantile regression for recurrent episode lengths\n")
  cat("  Quantiles: ", paste(x$tau, collapse = ", "), "\n", sep = "")
  cat("  Complete episodes assumed: ", if (isTRUE(x$complete)) "TRUE" else "FALSE", "\n", sep = "")
  if (!isTRUE(x$complete)) {
    cat("  Incomplete episodes dropped: ", x$n_dropped_incomplete, "\n", sep = "")
  }
  cat("  Time-dependent covariates: ",
      if (length(x$td_vars) == 0) "<none>" else paste(x$td_vars, collapse = ", "),
      "\n", sep = "")
  cat("  Time-independent covariates: ",
      if (length(x$z_vars) == 0) "<intercept only>" else paste(x$z_vars, collapse = ", "),
      "\n", sep = "")
  cat("  Episodes used: ", x$n_episodes_used, "\n", sep = "")
  cat("  Subjects observed (mi>0): ", x$n_subjects_observed, "\n", sep = "")
  if (!is.null(x$bootstrap) && !is.null(x$bootstrap$B)) {
    cat("  Bootstrap replicates: ", x$bootstrap$B, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.rqrecur <- function(object, ...) {
  co <- coef(object)
  out <- list(
    call = object$call,
    tau = object$tau,
    coef = co,
    complete = object$complete,
    n_dropped_incomplete = object$n_dropped_incomplete,
    td_vars = object$td_vars,
    z_vars = object$z_vars,
    n_episodes_used = object$n_episodes_used,
    n_subjects_observed = object$n_subjects_observed,
    censoring = object$censoring,
    bootstrap = object$bootstrap
  )
  class(out) <- "summary.rqrecur"
  out
}

#' @export
print.summary.rqrecur <- function(x, ...) {
  cat("<summary.rqrecur>\n")
  cat("Call:\n")
  print(x$call)
  cat("\nQuantiles: ", paste(x$tau, collapse = ", "), "\n", sep = "")
  cat("Complete episodes assumed: ", if (isTRUE(x$complete)) "TRUE" else "FALSE", "\n", sep = "")
  if (!isTRUE(x$complete)) {
    cat("Incomplete episodes dropped: ", x$n_dropped_incomplete, "\n", sep = "")
  }
  cat("Episodes used: ", x$n_episodes_used, "\n", sep = "")
  cat("Subjects observed (mi>0): ", x$n_subjects_observed, "\n", sep = "")
  cat("Censoring handling: ", x$censoring$mode, "\n\n", sep = "")
  if (!is.null(x$bootstrap) && !is.null(x$bootstrap$B)) {
    cat("Bootstrap replicates: ", x$bootstrap$B, "\n\n", sep = "")
  }
  cat("Coefficients (rows) by tau (columns):\n")
  print(round(x$coef, 6))
  invisible(x)
}

#' Extract coefficients
#'
#' @param object A fitted `"rqrecur"` object.
#' @param ... Unused.
#'
#' @return A numeric matrix with one column per `tau`.
#' @export
coef.rqrecur <- function(object, ...) {
  b1 <- object$coef$beta1
  b2 <- object$coef$beta2

  if (is.null(dim(b1))) {
    b1 <- matrix(b1, ncol = 1, dimnames = list(names(b1), as.character(object$tau[1])))
  }
  if (length(object$td_vars) == 0) {
    return(b1)
  }

  overlap <- intersect(rownames(b1), rownames(b2))
  if (length(overlap) > 0) {
    stop("Overlapping coefficient names between baseline and time-dependent parts: ",
         paste(overlap, collapse = ", "),
         "\nThis usually indicates the same variable was treated as both time-dependent and time-independent.",
         call. = FALSE)
  }
  rbind(b1, b2)
}

#' Predict episode-length quantiles for new covariates
#'
#' @param object A fitted `"rqrecur"` object.
#' @param newdata A data.frame containing the covariates used in the fit.
#' @param tau Optional subset of quantiles to predict. Defaults to all fitted `object$tau`.
#' @param ... Unused.
#'
#' @return A numeric matrix with one row per row of `newdata` and one column per `tau`.
#' @export
predict.rqrecur <- function(object, newdata, tau = NULL, ...) {
  if (missing(newdata) || is.null(newdata)) stop("`newdata` must be provided.", call. = FALSE)
  if (!is.data.frame(newdata)) stop("`newdata` must be a data.frame.", call. = FALSE)

  taus <- object$tau
  if (!is.null(tau)) {
    tau <- as.numeric(tau)
    keep <- tau %in% taus
    if (!all(keep)) stop("Some requested `tau` were not fitted: ", paste(tau[!keep], collapse = ", "), call. = FALSE)
    taus <- tau
  }
  tau_chr <- as.character(taus)

  # Baseline design matrix Z: reuse terms from the step-2 fit
  fit2_ref <- object$fit_step2[[1]]
  trm <- stats::delete.response(stats::terms(fit2_ref))
  mf <- stats::model.frame(trm, data = newdata, na.action = na.pass,
                           xlev = fit2_ref$xlevels)
  Z <- stats::model.matrix(trm, mf)

  # Time-dependent matrix L (numeric only)
  if (length(object$td_vars) > 0) {
    missing_td <- setdiff(object$td_vars, names(newdata))
    if (length(missing_td) > 0) stop("Missing time-dependent column(s) in newdata: ",
                                     paste(missing_td, collapse = ", "), call. = FALSE)
    L <- as.matrix(newdata[, object$td_vars, drop = FALSE])
    for (v in object$td_vars) {
      if (!is.numeric(newdata[[v]])) stop("Time-dependent covariate `", v, "` in newdata must be numeric.", call. = FALSE)
    }
    colnames(L) <- object$td_vars
  } else {
    L <- matrix(0, nrow(newdata), 0)
  }

  b1_mat <- object$coef$beta1[, tau_chr, drop = FALSE]
  b2_mat <- object$coef$beta2[, tau_chr, drop = FALSE]

  pred <- matrix(NA_real_, nrow = nrow(newdata), ncol = length(taus),
                 dimnames = list(NULL, tau_chr))

  for (k in seq_along(taus)) {
    b1 <- b1_mat[, k]
    names(b1) <- rownames(b1_mat)
    if (!all(colnames(Z) %in% names(b1))) {
      missing <- setdiff(colnames(Z), names(b1))
      stop("Internal mismatch in baseline design matrix. Missing coefficient(s): ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    mu1 <- as.numeric(Z %*% b1[colnames(Z)])

    if (ncol(L) > 0) {
      b2 <- b2_mat[, k]
      names(b2) <- rownames(b2_mat)
      mu2 <- as.numeric(L %*% b2[colnames(L)])
    } else {
      mu2 <- 0
    }
    pred[, k] <- mu1 + mu2
  }

  pred
}
