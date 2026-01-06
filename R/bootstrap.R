# Bootstrap and inference utilities for recurQR

#' Bootstrap a fitted recurrent-episode quantile regression model
#'
#' This function performs a **cluster (subject-level)** bootstrap: subjects are resampled with
#' replacement and their episode rows are carried along.
#'
#' The bootstrap is used to form **pointwise** confidence intervals for the coefficient
#' functions \eqn{\beta(\tau)} across the fitted \eqn{\tau} grid.
#'
#' @param object A fitted `"rqrecur"` object.
#' @param data The original long-format data used for fitting.
#' @param B Number of bootstrap replicates.
#' @param seed Optional integer seed for reproducibility.
#' @param progress Logical; if `TRUE`, show a progress bar.
#' @param max_tries Maximum number of attempted refits (to recover from occasional bootstrap
#'   failures). Default is `5 * B`.
#' @param ... Unused.
#'
#' @return The input `object`, augmented with a `bootstrap` element.
#' @export
rqrecur_boot <- function(object,
                        data,
                        B = 200,
                        seed = NULL,
                        progress = TRUE,
                        max_tries = NULL,
                        ...) {

  if (!inherits(object, "rqrecur")) stop("`object` must be of class 'rqrecur'.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)

  B <- as.integer(B)
  if (!is.finite(B) || B <= 0) stop("`B` must be a positive integer.", call. = FALSE)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed)) {
      stop("`seed` must be a single finite number (or NULL).", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  id_col <- object$id
  if (!id_col %in% names(data)) stop("Column `", id_col, "` not found in `data`.", call. = FALSE)

  # If the fitted object required a status column, we need it to exist for bootstrap refits.
  if (isFALSE(object$complete)) {
    if (is.null(object$status) || !object$status %in% names(data)) {
      stop(
        "The model was fitted with `complete = FALSE`, but `object$status` is missing or not a column in `data`.\n",
        "Please refit with `status` as a column name in `data` (recommended), or refit using `boot = TRUE` inside rqrecur().",
        call. = FALSE
      )
    }
  }

  ids <- unique(as.character(data[[id_col]]))
  n_ids <- length(ids)
  if (n_ids == 0) stop("No subjects found in `data`.", call. = FALSE)

  ref_coef <- coef(object)
  ref_names <- rownames(ref_coef)
  tau_chr <- colnames(ref_coef)

  if (is.null(max_tries)) max_tries <- as.integer(5 * B)
  max_tries <- as.integer(max_tries)
  if (!is.finite(max_tries) || max_tries < B) max_tries <- as.integer(5 * B)

  boot_array <- array(NA_real_, dim = c(length(ref_names), length(tau_chr), B),
                      dimnames = list(ref_names, tau_chr, paste0("b", seq_len(B))))

  pb <- NULL
  if (isTRUE(progress)) {
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit({
      try(utils::close(pb), silent = TRUE)
    }, add = TRUE)
  }

  b <- 0
  tries <- 0
  while (b < B && tries < max_tries) {
    tries <- tries + 1

    # Sample subjects with replacement
    draw <- sample(ids, size = n_ids, replace = TRUE)

    # Bind their rows and relabel IDs so duplicates remain separate clusters
    pieces <- vector("list", length(draw))
    for (k in seq_along(draw)) {
      d0 <- data[data[[id_col]] == draw[k], , drop = FALSE]
      d0[[id_col]] <- paste0(draw[k], "__boot", k)
      pieces[[k]] <- d0
    }
    dat_b <- do.call(rbind, pieces)

    # Stabilize factor levels if they existed in the original fit
    dat_b <- .apply_xlevels(dat_b, object$xlevels)

    fit_b <- tryCatch(
      rqrecur(
        formula = object$formula,
        data = dat_b,
        id = object$id,
        tau = object$tau,
        td_vars = object$td_vars,
        episode = object$episode,
        complete = object$complete,
        status = object$status,
        status_complete = object$status_complete,
        censor_time = object$censoring$censor_time_arg,
        w1 = object$w1_arg,
        eps = object$eps,
        boot = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(fit_b)) next

    co_b <- coef(fit_b)
    mat <- matrix(NA_real_, nrow = length(ref_names), ncol = length(tau_chr),
                  dimnames = list(ref_names, tau_chr))

    common <- intersect(rownames(co_b), ref_names)
    common_tau <- intersect(colnames(co_b), tau_chr)
    if (length(common) == 0 || length(common_tau) == 0) next

    mat[common, common_tau] <- co_b[common, common_tau, drop = FALSE]

    b <- b + 1
    boot_array[, , b] <- mat

    if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
  }

  if (b < B) {
    warning(
      "Only ", b, " successful bootstrap refits out of requested B=", B, ". ",
      "(Tried ", tries, " times.)\n",
      "This can happen if a bootstrap sample contains too few subjects with >=2 complete episodes ",
      "to estimate time-dependent effects.",
      call. = FALSE
    )
    boot_array <- boot_array[, , seq_len(b), drop = FALSE]
  }

  object$bootstrap <- list(
    B = dim(boot_array)[3],
    seed = seed,
    tries = tries,
    coef = boot_array
  )

  object
}

#' Pointwise confidence intervals for coefficient functions
#'
#' Computes pointwise bootstrap confidence intervals for \eqn{\beta(\tau)}.
#'
#' @param object A fitted `"rqrecur"` object with bootstrap results (see [rqrecur_boot()] or
#'   fit with `boot = TRUE` in [rqrecur()]).
#' @param parm Optional coefficient name(s). If `NULL`, returns all coefficients.
#' @param level Confidence level (default 0.95).
#' @param ... Unused.
#'
#' @return A list with elements `estimate`, `lower`, `upper`, and `level`.
#' @export
confint.rqrecur <- function(object, parm = NULL, level = 0.95, ...) {
  if (!inherits(object, "rqrecur")) stop("`object` must be of class 'rqrecur'.", call. = FALSE)
  if (is.null(object$bootstrap) || is.null(object$bootstrap$coef)) {
    stop("No bootstrap results found. Fit with `boot = TRUE` in `rqrecur()` (or call `rqrecur_boot()` first).", call. = FALSE)
  }

  level <- as.numeric(level)
  if (!is.finite(level) || level <= 0 || level >= 1) stop("`level` must be in (0, 1).", call. = FALSE)
  alpha <- (1 - level) / 2

  est <- coef(object)
  boot <- object$bootstrap$coef

  if (!is.null(parm)) {
    parm <- as.character(parm)
    missing <- setdiff(parm, rownames(est))
    if (length(missing) > 0) stop("Unknown coefficient name(s): ", paste(missing, collapse = ", "), call. = FALSE)
    est <- est[parm, , drop = FALSE]
    boot <- boot[parm, , , drop = FALSE]
  }

  # boot is [coef, tau, B]
  lower <- apply(boot, c(1, 2), stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE)
  upper <- apply(boot, c(1, 2), stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)

  list(estimate = est, lower = lower, upper = upper, level = level)
}
