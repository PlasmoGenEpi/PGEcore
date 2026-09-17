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
#' Collapses a per-parameter convergence table (one row per sampled
#' parameter, with R-hat and effective sample size) into a single row of
#' run-level health indicators. The per-parameter table answers "did this
#' parameter mix?"; this summary answers "did this run mix?" so runs can be
#' compared across populations, parameter sets, and tools without reading
#' hundreds of rows.
#'
#' @section Why one row per run:
#' Two failure modes matter in practice and are hard to see per parameter.
#' A chain that never moves off its initial state produces constant draws,
#' for which R-hat and ESS are undefined rather than bad. In a per-parameter
#' table that shows up as scattered `NA`s and is easy to miss. In this
#' summary it shows up as a high `pct_degenerate`, which is the standing
#' guard for frozen chains across the parallel-chain wrappers. The second
#' failure is broad, mild non-convergence, where no single parameter is
#' alarming but a large fraction sit above R-hat 1.01. The `pct_rhat_*`
#' columns expose that as a proportion rather than a worst case.
#'
#' @section How to read the columns:
#' * **Degenerate parameters.** A parameter is degenerate when its draws do
#'   not vary (`sd` of zero) or its R-hat is not finite. This can be benign,
#'   for example a COI fixed at 1 for a monoclonal specimen, or it can mean
#'   the sampler is stuck. Interpret `pct_degenerate` against what the model
#'   is expected to hold constant. A value near 100% means the chain did not
#'   move, whatever the R-hat columns say.
#' * **R-hat.** Computed on usable parameters only. `max_rhat` is the worst
#'   case. The three `pct_rhat_gt_*` columns give the share of usable
#'   parameters above 1.01 (the strict modern threshold), 1.05 (a common
#'   working tolerance), and 1.1 (the older lenient cutoff). A well-mixed run
#'   has a `max_rhat` near 1 and small percentages in all three.
#' * **Effective sample size.** `ess_bulk` reflects how well posterior means
#'   and medians are estimated. `ess_tail` reflects the quantiles that feed
#'   credible intervals. The minima flag the worst parameter; the median
#'   describes the typical one. A rule of thumb is at least 100 effective
#'   draws per chain for a stable estimate. ESS can be `NA` for a parameter
#'   whose R-hat is finite, and those are skipped rather than collapsing the
#'   column to `NA`.
#'
#' Because R-hat and ESS are computed on usable parameters only, always read
#' them alongside `pct_degenerate`. A run with many degenerate parameters can
#' look clean on `max_rhat` simply because the parameters that would have
#' failed were excluded. `pct_*` values are percentages of `n_usable`, or
#' `NA` when nothing is usable.
#'
#' Convergence is a property of the sampler, not of the estimate. A run can
#' mix well and still be inaccurate, and the best-mixing parameter set is
#' not necessarily the most accurate one, so report these diagnostics as a
#' separate axis from accuracy rather than using them to pick a model.
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
