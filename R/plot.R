# Plotting utilities for recurQR

#' Plot fitted coefficient functions with optional pointwise confidence intervals
#'
#' This is a convenience method for plotting the estimated coefficient functions
#' \eqn{\beta(\tau)} across the fitted \eqn{\tau} grid.
#'
#' Confidence intervals are obtained from subject-level bootstrap samples (see
#' [rqrecur_boot()] or fit with `boot = TRUE` in [rqrecur()]).
#'
#' @param x A fitted `rqrecur` object.
#' @param parm Coefficient name(s) to plot. If `NULL`, plots all coefficients.
#' @param level Confidence level for pointwise intervals (default 0.95).
#' @param ci Logical; whether to draw confidence intervals when available.
#' @param ncol Number of columns in the multi-panel layout when plotting multiple
#'   coefficients. If `NULL`, a square-ish layout is chosen automatically.
#' @param ask Logical; if `TRUE`, pause between pages when many panels are drawn.
#' @param ... Passed to the underlying `plot()` call (e.g., `ylim`, `lwd`).
#'
#' @return Coefficients plot.
#' @export
plot.rqrecur <- function(x,
                         parm = NULL,
                         level = 0.95,
                         ci = TRUE,
                         ncol = NULL,
                         ask = NULL,
                         ...) {

  if (!inherits(x, "rqrecur")) stop("`x` must be of class 'rqrecur'.", call. = FALSE)

  est <- coef(x)
  taus <- as.numeric(colnames(est))

  if (is.null(parm)) {
    parm <- rownames(est)
  } else {
    parm <- as.character(parm)
    missing <- setdiff(parm, rownames(est))
    if (length(missing) > 0) stop("Unknown coefficient name(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  # Decide whether we can draw CIs
  ci_obj <- NULL
  if (isTRUE(ci)) {
    if (!is.null(x$bootstrap) && !is.null(x$bootstrap$coef)) {
      ci_obj <- confint(x, parm = parm, level = level)
    } else {
      warning("No bootstrap results found; plotting without confidence intervals.", call. = FALSE)
    }
  }

  n_plot <- length(parm)
  if (n_plot == 0) return(invisible(x))

  # Layout
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(n_plot))
  }
  nrow <- ceiling(n_plot / ncol)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (!is.null(ask)) {
    graphics::par(ask = isTRUE(ask))
  }
  graphics::par(mfrow = c(nrow, ncol))

  for (nm in parm) {
    y <- as.numeric(est[nm, ])

    if (!is.null(ci_obj)) {
      lo <- as.numeric(ci_obj$lower[nm, ])
      hi <- as.numeric(ci_obj$upper[nm, ])
      ylim <- range(c(y, lo, hi), finite = TRUE)
      graphics::plot(taus, y,
                     type = "l",
                     xlab = "Quantile level",
                     ylab = "Coefficient",
                     main = nm,
                     ylim = ylim,
                     ...)
      graphics::lines(taus, lo, lty = 3)
      graphics::lines(taus, hi, lty = 3)
    } else {
      graphics::plot(taus, y,
                     type = "l",
                     xlab = "Quantile level",
                     ylab = "Coefficient",
                     main = nm,
                     ...)
    }

    graphics::abline(h = 0, lty = 2, col = "gray")
  }

  invisible(x)
}
