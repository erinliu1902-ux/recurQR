# Internal helper utilities for recurQR

.check_scalar_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || nchar(x) == 0) {
    stop("`", name, "` must be a single non-empty string.", call. = FALSE)
  }
}

.check_cols_exist <- function(data, cols, where = "data") {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Missing column(s) in ", where, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

## Filter to complete episodes
##
## In the paper, m_i is defined as the number of *completely observed* episodes.
## Users may provide a long dataset that also contains an incomplete (censored)
## final episode. The `complete` flag controls whether we need to filter those
## rows out.
##
## If `complete == TRUE`, all rows are assumed complete and kept.
## If `complete == FALSE`, we keep only rows where `status` indicates a complete
## episode.
.filter_complete_episodes <- function(data,
                                     complete,
                                     status = NULL,
                                     status_complete = 1) {
  if (!is.logical(complete) || length(complete) != 1 || is.na(complete)) {
    stop("`complete` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  n <- nrow(data)
  if (n == 0) return(logical(0))
  if (complete) return(rep(TRUE, n))

  # Need a censoring/episode-status indicator.
  if (is.null(status)) {
    # Best-effort auto-detect
    cand <- intersect(c("status", "delta", "event", "episode_complete", "complete", "censored"), names(data))
    if (length(cand) == 0) {
      stop(
        "`complete = FALSE` requires an episode censoring indicator via `status` (column name or vector).\n",
        "Typical choices: status=\"delta\" (1=complete,0=censored) or status=\"status\".",
        call. = FALSE
      )
    }
    status <- cand[1]
  }

  if (is.character(status) && length(status) == 1) {
    if (!status %in% names(data)) stop("Column `", status, "` not found in `data`.", call. = FALSE)
    v <- data[[status]]
  } else if (is.logical(status) && length(status) == n) {
    v <- status
  } else if (is.numeric(status) && length(status) == n) {
    v <- status
  } else {
    stop("`status` must be a single column name or a vector of length nrow(data).", call. = FALSE)
  }

  keep <- NULL
  if (is.logical(v)) {
    keep <- v
  } else if (is.numeric(v)) {
    keep <- as.numeric(v) == status_complete
  } else {
    stop("`status` must be logical or numeric (e.g., 1=complete, 0=censored).", call. = FALSE)
  }

  keep[is.na(keep)] <- FALSE
  keep
}

## Apply factor levels captured from a reference fit (helps bootstrap stability)
.apply_xlevels <- function(data, xlevels) {
  if (is.null(xlevels) || length(xlevels) == 0) return(data)
  for (nm in names(xlevels)) {
    if (nm %in% names(data)) {
      data[[nm]] <- factor(data[[nm]], levels = xlevels[[nm]])
    }
  }
  data
}

.resolve_episode_order <- function(data, id, episode = NULL) {
  if (!is.null(episode)) {
    .check_scalar_string(episode, "episode")
    if (!episode %in% names(data)) stop("Column `", episode, "` not found in `data`.", call. = FALSE)
    ep <- data[[episode]]
    if (!is.numeric(ep)) stop("`episode` column must be numeric/integer.", call. = FALSE)
    return(as.integer(ep))
  }

  if ("indx" %in% names(data) && is.numeric(data[["indx"]])) {
    return(as.integer(data[["indx"]]))
  }
  if ("episode" %in% names(data) && is.numeric(data[["episode"]])) {
    return(as.integer(data[["episode"]]))
  }

  # try to order by 'end' time within id (treating 'end' as start-of-episode time)
  if ("end" %in% names(data) && is.numeric(data[["end"]])) {
    # create within-id ranks of end (ties broken by row order)
    ord <- order(data[[id]], data[["end"]], seq_len(nrow(data)))
    rk <- integer(nrow(data))
    rk[ord] <- ave(data[["end"]][ord], as.character(data[[id]][ord]),
                   FUN = function(x) seq_along(x))
    return(as.integer(rk))
  }

  # fallback: row order within id
  ord <- order(data[[id]], seq_len(nrow(data)))
  rk <- integer(nrow(data))
  rk[ord] <- ave(seq_along(ord), as.character(data[[id]][ord]),
                 FUN = function(x) seq_along(x))
  as.integer(rk)
}

.resolve_w1_name <- function(data, w1 = NULL) {
  if (!is.null(w1)) {
    .check_scalar_string(w1, "w1")
    if (!w1 %in% names(data)) stop("Column `", w1, "` not found in `data`.", call. = FALSE)
    return(w1)
  }
  if ("T1" %in% names(data)) return("T1")
  if ("w1" %in% names(data)) return("w1")
  if ("Wi1" %in% names(data)) return("Wi1")
  if ("end" %in% names(data)) return(NULL) # will compute from end of first episode
  stop("Unable to infer W_i1. Please provide `w1` (column name) or include a `T1`/`w1` column.", call. = FALSE)
}

.compute_w1_per_row <- function(data, id, episode_col, w1_name = NULL) {
  id_chr <- as.character(data[[id]])
  first <- !duplicated(id_chr[order(id_chr, data[[episode_col]])])
  # reorder to get first per subject correctly
  ord <- order(id_chr, data[[episode_col]])
  first_ord <- !duplicated(id_chr[ord])
  idx_first <- ord[first_ord]

  if (!is.null(w1_name)) {
    if (!is.numeric(data[[w1_name]])) stop("`", w1_name, "` must be numeric.", call. = FALSE)
    w1_by_id <- data[[w1_name]][idx_first]
  } else {
    # infer as 'end' time of first episode
    if (!("end" %in% names(data) && is.numeric(data[["end"]]))) {
      stop("Cannot infer W_i1 without an `end` column.", call. = FALSE)
    }
    w1_by_id <- data[["end"]][idx_first]
  }
  names(w1_by_id) <- id_chr[idx_first]
  w1_by_id[id_chr]
}

# Empirical survival function G_hat(t) = P(C >= t)
.emp_survival <- function(c_times, t) {
  c_sorted <- sort(as.numeric(c_times))
  n <- length(c_sorted)
  if (n == 0) stop("Empty censoring times.", call. = FALSE)

  t <- as.numeric(t)
  # Count of C < t (strict), matching the paper's empirical survival
  #   S_C(t) = #{ C >= t } / n
  # exactly (including ties). Using `left.open = TRUE` makes findInterval
  # return the number of elements strictly less than `t`.
  n_lt <- findInterval(t, c_sorted, left.open = TRUE)
  (n - n_lt) / n
}
