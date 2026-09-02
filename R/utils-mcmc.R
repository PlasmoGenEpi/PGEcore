#' Summarize posterior draws with convergence diagnostics
#'
#' Keeps the subset of [posterior::summarise_draws()] columns reported by the
#' MCMC wrappers, so convergence tables share one schema across tools.
#'
#' @param draws A `draws` object (for example from
#'   [posterior::as_draws_array()]).
#' @return A data frame with one row per parameter and the columns `variable`,
#'   `mean`, `median`, `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, and `ess_tail`.
#' @keywords internal
summarize_convergence_draws <- function(draws) {
  check_suggested_pkg("posterior", "MCMC convergence diagnostics")
  summary_tbl <- posterior::summarise_draws(draws)
  keep <- c(
    "variable", "mean", "median", "sd", "q5", "q95",
    "rhat", "ess_bulk", "ess_tail"
  )
  as.data.frame(summary_tbl[, keep])
}

#' Reduce a convergence table to one run-level row
#'
#' A parameter is degenerate when its draws do not vary (`sd` of zero) or its
#' R-hat is not finite; R-hat and ESS are undefined for those and they are
#' excluded from the R-hat and ESS columns below. `pct_*` values are percentages
#' of `n_usable`, or `NA` when nothing is usable.
#'
#' @param convergence A convergence table from
#'   [summarize_convergence_draws()].
#' @param group_cols Optional column names to summarize within, for tables that
#'   stack several fits (for example FEM's `group_id`).
#' @return A data frame with the columns `n_variables`, `n_degenerate`,
#'   `pct_degenerate`, `n_usable`, `max_rhat`, `pct_rhat_gt_1_01`,
#'   `pct_rhat_gt_1_05`, `pct_rhat_gt_1_1`, `min_ess_bulk`, `median_ess_bulk`
#'   and `min_ess_tail`, preceded by any `group_cols`.
#' @keywords internal
summarize_convergence_run <- function(convergence, group_cols = NULL) {
  summarize_one <- function(df) {
    degenerate <- is.na(df$sd) | df$sd == 0 |
      is.na(df$rhat) | !is.finite(df$rhat)
    ok <- df[!degenerate, , drop = FALSE]
    n_usable <- nrow(ok)
    pct <- function(x) {
      if (n_usable == 0L) NA_real_ else 100 * mean(x)
    }
    # ESS can be NA for a parameter whose R-hat is finite, so those are skipped
    # rather than collapsing the whole column to NA.
    min_or_na <- function(x) {
      if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
    }
    data.frame(
      n_variables = nrow(df),
      n_degenerate = sum(degenerate),
      pct_degenerate = if (nrow(df) == 0L) NA_real_ else 100 * mean(degenerate),
      n_usable = n_usable,
      max_rhat = if (n_usable == 0L) NA_real_ else max(ok$rhat),
      pct_rhat_gt_1_01 = pct(ok$rhat > 1.01),
      pct_rhat_gt_1_05 = pct(ok$rhat > 1.05),
      pct_rhat_gt_1_1 = pct(ok$rhat > 1.1),
      min_ess_bulk = min_or_na(ok$ess_bulk),
      median_ess_bulk = if (all(is.na(ok$ess_bulk))) {
        NA_real_
      } else {
        stats::median(ok$ess_bulk, na.rm = TRUE)
      },
      min_ess_tail = min_or_na(ok$ess_tail)
    )
  }

  if (is.null(group_cols)) {
    return(summarize_one(convergence))
  }
  keys <- unique(convergence[, group_cols, drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    sel <- rep(TRUE, nrow(convergence))
    for (g in group_cols) {
      sel <- sel & convergence[[g]] == keys[[g]][i]
    }
    cbind(keys[i, , drop = FALSE], summarize_one(convergence[sel, , drop = FALSE]))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
