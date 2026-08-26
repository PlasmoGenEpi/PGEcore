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
