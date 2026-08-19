suppressPackageStartupMessages(library(dplyr))

#' Convert prev df to STAVE Format for variant column
#'
#' This function processes a data frame by generating a `variant` column, 
#' which combines `gene_id`, `aa_position`, and `aa` into a single identifier.
#' It then returns only the `variant` and columns specified using `additional_columns``.
#'
#' @param df A data frame containing at least the columns:
#'   - `gene_id` (character): Gene identifier
#'   - `aa_position` (numeric or integer): Amino acid position
#'   - `aa` (character): Amino acid at the given position
#'   - any `additional_columns` specified
#'
#' @return A data frame with columns:
#'   - `variant` (character): Concatenated string of `gene_id:aa_position:aa`
#'   - any `additional_columns` specified
#'
#' @examples
#' df <- data.frame(
#'   gene_id = c("PF3D7_0417200.1", "PF3D7_0417200.1"),
#'   aa_position = c(51, 59),
#'   aa = c("I", "R"),
#'   prev = c(0.5, 1)
#' )
#' convert_to_stave(df)
convert_single_locus_table_to_stave <- function(df, additional_columns=NULL) {

  df %>%
    ungroup %>%
    mutate(variant = paste(gene_id, aa_position, aa, sep = ":")) %>%
    {
      if (is.null(additional_columns)) {
        select(., variant)
      } else {
        select(., variant, all_of(additional_columns))
      }
    }

}

#' Summarize MCMC convergence diagnostics from a posterior draws array
#'
#' Wraps `posterior::summarise_draws()` to produce the convergence diagnostics
#' table used by the MCMC wrapper scripts. Parameters that never move across
#' the retained samples yield `NaN` for `rhat`/`ess_bulk`/`ess_tail`; this means
#' "did not move", not a failure, and is passed through unchanged.
#'
#' @param draws A `posterior::draws_array` (or any object accepted by
#'   `posterior::summarise_draws()`) with one variable per estimated parameter.
#'
#' @return A data frame with one row per parameter and the columns: variable,
#'   mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail.
summarize_convergence_draws <- function(draws) {
  summary_tbl <- posterior::summarise_draws(draws)
  keep <- c(
    "variable", "mean", "median", "sd", "q5", "q95",
    "rhat", "ess_bulk", "ess_tail"
  )
  as.data.frame(summary_tbl[, keep])
}
