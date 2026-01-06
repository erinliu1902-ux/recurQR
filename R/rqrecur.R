#' Fit the quantile regression model for recurrent episode lengths
#'
#' Implements the two-step estimation procedure described in
#' *Exploring the Heterogeneity in Recurrent Episode Lengths Based On Quantile Regression*
#' (Liu, Umpierrez, and Peng).
#'
#' The method targets the conditional quantiles of episode lengths:
#'
#' \deqn{Q_{X_{ij}}(\tau \mid Z_i, L_{ij}) = \beta_1(\tau)^\top Z_i + \beta_2(\tau)^\top L_{ij}}
#'
#' and addresses dependent truncation, dependent censoring, and informative cluster size
#' via a two-step strategy:
#' 1) estimate \eqn{\beta_2(\tau)} from within-subject differences, and
#' 2) estimate \eqn{\beta_1(\tau)} by weighted quantile regression.
#'
#' @param formula A formula of the form `episode_length ~ x1 + x2 + ...` including both
#'   time-independent and time-dependent covariates. Currently only additive terms are supported
#'   (no interactions or transformations).
#' @param data A data.frame in long format with one row per (observed) episode.
#' @param id Name of the subject id column.
#' @param tau Quantile level(s) in `(0, 1)`. A numeric vector.
#' @param td_vars Character vector giving the *time-dependent* covariate names among the predictors
#'   in `formula`. The remaining predictors are treated as time-independent.
#'   Use `character(0)` if there are no time-dependent covariates.
#' @param episode Optional name of an episode-order column. If `NULL`, the function will use row order within subject.
#' @param complete Logical. If `TRUE`, `data` is assumed to contain only complete episodes.
#'   If `FALSE`, incomplete (censored) episodes are removed using `status`.
#' @param status Episode censoring indicator used when `complete = FALSE`. Either a column name
#'   in `data` or a logical/numeric vector of length `nrow(data)`. For a numeric
#'   indicator, `status_complete` specifies which value corresponds to a complete episode.
#' @param status_complete Value in `status` that indicates a complete episode when `status` is
#'   numeric. Default is `1` (1=complete, 0=censored).
#' @param censor_time Censoring/follow-up times.
#'   - If a numeric vector, it is treated as the cohort censoring times \eqn{C_1,\dots,C_n}.
#'   - If a string, it is interpreted as a column in `data`; the function will extract it from `data`.
#'   If `NULL`, the function will look for a `censor_time` column in `data`. 
#' @param w1 Name of the column containing \eqn{W_{i1}} (time from study entry to start of
#'   the first episode). If `NULL`, the function tries `"T1"` then `"w1"`, then uses the `"end"`
#'   time of the first episode if available. If no truncation correction is performed, `w1` is not required.
#' @param eps Small positive constant used to avoid division by zero in `1/G(t)`. Default `1e-8`.
#' @param boot Logical. If `TRUE`, run a **subject-level** bootstrap inside `rqrecur()` and attach
#'   pointwise confidence intervals to the returned object.
#' @param B Number of bootstrap replicates (only used when `boot = TRUE`). Default is `200`.
#' @param bootstrap Deprecated alias for `B`. If provided, it overrides `B` and also sets
#'   `boot = (bootstrap > 0)` for backward compatibility.
#' @param seed Optional integer seed for reproducible bootstrapping.
#' @param progress Logical; if `TRUE`, show a progress bar during bootstrapping.
#' @param ... Currently unused (reserved).
#'
#' @return An object of class `"rqrecur"` containing coefficient estimates for each `tau`.
#'
#' @examples
#' dat <- read_example_sim_data()
#' fit <- rqrecur(
#'   epi_length ~ z1 + z2 + z3,
#'   data = dat,
#'   id = "id",
#'   td_vars = "z3",
#'   tau = c(0.25, 0.5, 0.75),
#'   # In real use, provide censor_time from your full cohort
#'   censor_time = NULL
#' )
#' # Predict quantiles for new covariates
#' newdat <- data.frame(z1 = 1.3, z2 = 1, z3 = 5)
#' predict(fit, newdata = newdat)
#'
#' @export
rqrecur <- function(formula,
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
                    ...) {

  .check_scalar_string(id, "id")
  if (!inherits(formula, "formula")) stop("`formula` must be a formula.", call. = FALSE)
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)

  # --- backward-compatible alias: `censor_times` (deprecated)
  extra_args <- list(...)
  if (!is.null(extra_args$censor_times)) {
    warning("`censor_times` is deprecated. Please use `censor_time` instead.", call. = FALSE)
    if (is.null(censor_time)) {
      censor_time <- extra_args$censor_times
    }
  }

  tau <- as.numeric(tau)
  if (any(!is.finite(tau)) || any(tau <= 0 | tau >= 1)) {
    stop("`tau` must be in (0, 1).", call. = FALSE)
  }
  if (length(unique(tau)) != length(tau)) {
    tau <- sort(unique(tau))
  }

  if (missing(td_vars)) stop("`td_vars` must be provided (possibly character(0)).", call. = FALSE)
  if (!is.character(td_vars)) stop("`td_vars` must be a character vector.", call. = FALSE)

  # --- validate formula terms (we keep this strict to make 'delta' unambiguous)
  trm <- terms(formula)
  if (attr(trm, "intercept") != 1) {
    stop("The model must include an intercept (do not use -1 / 0 + in the formula).", call. = FALSE)
  }
  term_labels <- attr(trm, "term.labels")
  bad <- term_labels[!grepl("^[A-Za-z.][A-Za-z0-9._]*$", term_labels)]
  if (length(bad) > 0) {
    stop(
      "Only simple additive variable terms are supported in `formula` for now.\n",
      "Problematic term(s): ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  # Ensure td_vars are in the formula
  if (length(td_vars) > 0) {
    missing_td <- setdiff(td_vars, term_labels)
    if (length(missing_td) > 0) {
      stop("The following `td_vars` are not present in `formula`: ",
           paste(missing_td, collapse = ", "), call. = FALSE)
    }
  }
  z_vars <- setdiff(term_labels, td_vars)

  # --- handle completeness (drop incomplete/censored episodes if needed)
  dat <- data
  if (!id %in% names(dat)) stop("Column `", id, "` not found in `data`.", call. = FALSE)

  keep <- .filter_complete_episodes(dat,
                                   complete = complete,
                                   status = status,
                                   status_complete = status_complete)
  n_dropped_incomplete <- sum(!keep)
  dat <- dat[keep, , drop = FALSE]
  if (nrow(dat) == 0) stop("No complete episodes remain after filtering.", call. = FALSE)

  # --- resolve episode order
  ep <- .resolve_episode_order(dat, id = id, episode = episode)
  dat[[".rqrecur_episode"]] <- ep

  # --- extract response name and values
  # we use model.frame to respect NA handling consistent with the formula
  mf_full <- model.frame(formula, data = dat, na.action = na.omit)
  # after na.omit, we need to subset `dat` the same way
  dat <- dat[rownames(mf_full), , drop = FALSE]
  y <- model.response(mf_full)
  y_name <- as.character(formula[[2]])

  if (!is.numeric(y)) stop("Response must be numeric episode length.", call. = FALSE)

  # --- resolve censoring times (vector for G-hat)
  G_mode <- "none"
  c_times <- NULL

  # If `censor_time` is NULL, fall back to a `censor_time` column in the data (if present).
  #
  # Note: in long-format episode data, `censor_time` is typically repeated on every episode row
  # (dimension sum_i m_i). We take one value per subject when estimating G(t).
  if (is.null(censor_time)) {
    if ("censor_time" %in% names(dat)) {
      censor_time <- "censor_time"
    } else {
      warning(
        "`censor_time` not provided and no `censor_time` column found in `data`. ",
        "Setting G(t) = 1 (no truncation correction).",
        call. = FALSE
      )
      G_mode <- "none"
    }
  }

  if (is.null(censor_time)) {
    # no truncation correction
    G_mode <- "none"
  } else if (is.character(censor_time) && length(censor_time) == 1) {
    if (!censor_time %in% names(dat)) stop("Column `", censor_time, "` not found in `data`.", call. = FALSE)

    # One censor time per subject (best effort). Warn if the column is not constant within a subject.
    tmp <- dat[, c(id, censor_time)]
    tmp <- tmp[!is.na(tmp[[censor_time]]), , drop = FALSE]

    multi <- tapply(tmp[[censor_time]], tmp[[id]], function(v) {
      v <- v[is.finite(v)]
      length(unique(v)) > 1
    })
    if (any(multi, na.rm = TRUE)) {
      warning(
        "Multiple distinct `", censor_time, "` values found within some subjects; ",
        "using the first observed value per subject.",
        call. = FALSE
      )
    }

    c_by_id <- tapply(tmp[[censor_time]], tmp[[id]], function(v) {
      v <- v[is.finite(v)]
      if (length(v) == 0) return(NA_real_)
      v[1]
    })

    c_times <- as.numeric(c_by_id)
    c_times <- c_times[is.finite(c_times)]
    if (length(c_times) == 0) stop("No finite censoring times found.", call. = FALSE)
    G_mode <- "from_column"
  } else if (is.numeric(censor_time)) {
    c_times <- as.numeric(censor_time)
    c_times <- c_times[is.finite(c_times)]
    if (length(c_times) == 0) stop("No finite `censor_time` values provided.", call. = FALSE)
    G_mode <- "from_vector"
  } else {
    stop("`censor_time` must be NULL, a numeric vector, or a single column name.", call. = FALSE)
  }

  # --- resolve W_i1 (only needed when truncation correction is used)
  w1_vec <- rep(0, nrow(dat))
  w1_name <- NULL
  if (G_mode != "none") {
    w1_name <- .resolve_w1_name(dat, w1 = w1)
    w1_vec <- .compute_w1_per_row(dat, id = id, episode_col = ".rqrecur_episode", w1_name = w1_name)
    if (any(!is.finite(w1_vec))) {
      stop("Non-finite W_i1 values encountered. Please check your `w1`/time columns.", call. = FALSE)
    }
  }

  # --- build per-subject summaries: first episode values, cluster size, L1
  ord2 <- order(dat[[id]], dat[[".rqrecur_episode"]])
  dat <- dat[ord2, , drop = FALSE]
  y <- y[ord2]
  w1_vec <- w1_vec[ord2]
  id_chr <- as.character(dat[[id]])

  first_row <- !duplicated(id_chr)
  ids_obs <- unique(id_chr)

  # IMPORTANT: do *not* overwrite the names coming from `table()`.
  # `id_chr` is a character vector; `table(id_chr)` is ordered lexicographically.
  # If ids are numeric (e.g., 1,2,...,10,11,...) then `unique(id_chr)` follows the
  # current row order (typically numeric order after sorting), which can differ
  # from the lexicographic order used by `table()`. Overwriting names would then
  # misalign cluster sizes across subjects and corrupt weights.
  tab_mi <- table(id_chr)
  mi_by_id <- as.integer(tab_mi)
  names(mi_by_id) <- names(tab_mi)

  # first episode length Xi1
  xi1_by_id <- y[first_row]
  names(xi1_by_id) <- id_chr[first_row]
  xi1_row <- xi1_by_id[id_chr]

  # L and L1 (time-dependent)
  if (length(td_vars) > 0) {
    .check_cols_exist(dat, td_vars, where = "`data`")
    # enforce numeric td vars (delta needs numeric)
    for (v in td_vars) {
      if (!is.numeric(dat[[v]])) stop("Time-dependent covariate `", v, "` must be numeric.", call. = FALSE)
    }
    L_mat <- as.matrix(dat[, td_vars, drop = FALSE])
    colnames(L_mat) <- td_vars
    L1_by_id <- L_mat[first_row, , drop = FALSE]
    rownames(L1_by_id) <- id_chr[first_row]
    L1_mat <- L1_by_id[id_chr, , drop = FALSE]
    dL_mat <- L_mat - L1_mat
  } else {
    L_mat <- matrix(0, nrow(dat), 0)
    L1_mat <- matrix(0, nrow(dat), 0)
    dL_mat <- matrix(0, nrow(dat), 0)
  }

  # --- Step 1: estimate beta2(tau) from within-subject differences
  step1_fit <- NULL
  beta2_hat <- matrix(0, nrow = length(td_vars), ncol = length(tau))
  rownames(beta2_hat) <- td_vars
  colnames(beta2_hat) <- as.character(tau)
  eta_hat <- rep(NA_real_, length(tau))
  names(eta_hat) <- as.character(tau)

  if (length(td_vars) > 0) {
    is_first <- first_row
    mi_row <- mi_by_id[id_chr]
    keep1 <- (!is_first) & (mi_row > 1)
    if (sum(keep1) == 0) {
      stop("No subjects with >=2 complete episodes; cannot estimate time-dependent effects (beta2).", call. = FALSE)
    }

    dX <- y - xi1_row
    d_df <- data.frame(dX = dX[keep1], stringsAsFactors = FALSE)
    for (k in seq_along(td_vars)) {
      d_df[[paste0("d_", td_vars[k])]] <- dL_mat[keep1, k]
    }
    f1_rhs <- paste0("d_", td_vars)
    f1 <- stats::as.formula(paste("dX ~", paste(f1_rhs, collapse = " + ")))
    step1_fit <- quantreg::rq(f1, data = d_df, tau = tau)

    c1 <- stats::coef(step1_fit)
    if (is.null(dim(c1))) {
      c1 <- matrix(c1, ncol = 1, dimnames = list(names(stats::coef(step1_fit)), as.character(tau[1])))
    }
    eta_hat <- c1["(Intercept)", , drop = TRUE]
    for (k in seq_along(td_vars)) {
      beta2_hat[k, ] <- c1[paste0("d_", td_vars[k]), ]
    }
  }

  # --- Step 2: estimate beta1(tau) by weighted QR of adjusted outcomes on Z
  # Build baseline design via a baseline-only formula so prediction can reuse it.
  if (length(z_vars) > 0) {
    .check_cols_exist(dat, z_vars, where = "`data`")
    f2_rhs <- paste(z_vars, collapse = " + ")
  } else {
    f2_rhs <- "1"
  }
  f2 <- stats::as.formula(paste("y_adj ~", f2_rhs))

  step2_fits <- vector("list", length(tau))
  names(step2_fits) <- as.character(tau)

  beta1_hat_list <- vector("list", length(tau))
  names(beta1_hat_list) <- as.character(tau)

  mi_row <- mi_by_id[id_chr]

  for (k_tau in seq_along(tau)) {
    b2 <- if (length(td_vars) > 0) beta2_hat[, k_tau] else numeric(0)
    # adjusted outcome: X - beta2' L
    if (length(td_vars) > 0) {
      y_adj <- as.numeric(y - (L_mat %*% b2))
    } else {
      y_adj <- y
    }

    # weights: 1/mi * 1/G( W1 + X - beta2' dL )
    if (G_mode == "none") {
      w <- 1 / mi_row
    } else {
      t_eval <- if (length(td_vars) > 0) {
        as.numeric(w1_vec + y - (dL_mat %*% b2))
      } else {
        as.numeric(w1_vec + y)
      }
      G_hat <- .emp_survival(c_times, t_eval)
      w <- (1 / mi_row) * (1 / pmax(G_hat, eps))
    }

    step2_df <- dat
    step2_df[["y_adj"]] <- y_adj
    step2_fit <- quantreg::rq(f2, data = step2_df, tau = tau[k_tau], weights = w)
    step2_fits[[k_tau]] <- step2_fit
    beta1_hat_list[[k_tau]] <- stats::coef(step2_fit)
  }

  # combine beta1 into a matrix (rows = coef names)
  beta1_names <- unique(unlist(lapply(beta1_hat_list, names)))
  beta1_hat <- matrix(NA_real_, nrow = length(beta1_names), ncol = length(tau),
                      dimnames = list(beta1_names, as.character(tau)))
  for (k_tau in seq_along(tau)) {
    b <- beta1_hat_list[[k_tau]]
    beta1_hat[names(b), k_tau] <- b
  }

  status_col <- if (is.character(status) && length(status) == 1) status else NULL

  out <- list(
    call = match.call(),
    formula = formula,
    response = y_name,
    id = id,
    episode = episode,
    complete = complete,
    status = status_col,
    status_complete = status_complete,
    n_dropped_incomplete = n_dropped_incomplete,
    tau = tau,
    td_vars = td_vars,
    z_vars = z_vars,
    coef = list(beta1 = beta1_hat, beta2 = beta2_hat),
    eta = eta_hat,
    fit_step1 = step1_fit,
    fit_step2 = step2_fits,
    xlevels = step2_fits[[1]]$xlevels,
    censoring = list(
      mode = G_mode,
      censor_time_n = ifelse(is.null(c_times), NA_integer_, length(c_times)),
      censor_time_arg = censor_time
    ),
    w1_name = w1_name,
    w1_arg = w1,
    eps = eps,
    n_subjects_observed = length(ids_obs),
    n_episodes_used = nrow(dat)
  )
  class(out) <- "rqrecur"

  # Optional bootstrap for pointwise confidence intervals
  #
  # New interface (v0.1.2): `boot` (TRUE/FALSE) + `B`.
  # Backward compatible alias: `bootstrap`.
  if (!is.null(bootstrap)) {
    warning("`bootstrap` is deprecated. Please use `boot` and `B` instead.", call. = FALSE)
    bootstrap <- as.integer(bootstrap)
    if (!is.finite(bootstrap) || bootstrap < 0) {
      stop("`bootstrap` must be a nonnegative integer.", call. = FALSE)
    }
    B <- bootstrap
    boot <- (B > 0)
  }

  boot <- isTRUE(boot)
  if (!is.numeric(B) || length(B) != 1 || !is.finite(B)) {
    stop("`B` must be a single finite number (and will be coerced to integer).", call. = FALSE)
  }
  B <- as.integer(B)

  if (boot) {
    if (B <= 0) stop("`B` must be a positive integer when `boot = TRUE`.", call. = FALSE)
    out <- rqrecur_boot(out, data = data, B = B, seed = seed, progress = progress)
  }

  out
}