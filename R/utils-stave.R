#' Convert a single-locus table to STAVE-style variant identifiers
#'
#' Builds a `variant` column by concatenating `gene_id`, `aa_position`, and
#' `aa` as `gene_id:aa_position:aa`, then returns `variant` plus any requested
#' additional columns.
#'
#' @param df A data frame with columns `gene_id`, `aa_position`, and `aa`,
#'   plus any columns named in `additional_columns`.
#' @param additional_columns Optional character vector of extra columns to keep.
#'
#' @return A data frame with a `variant` column and any `additional_columns`.
#'
#' @examples
#' df <- data.frame(
#'   gene_id = c("PF3D7_0417200.1", "PF3D7_0417200.1"),
#'   aa_position = c(51, 59),
#'   aa = c("I", "R"),
#'   prev = c(0.5, 1)
#' )
#' convert_single_locus_table_to_stave(df, additional_columns = "prev")
#'
#' @export
convert_single_locus_table_to_stave <- function(df, additional_columns = NULL) {
  out <- df |>
    dplyr::ungroup() |>
    dplyr::mutate(
      variant = paste(.data$gene_id, .data$aa_position, .data$aa, sep = ":")
    )

  if (is.null(additional_columns)) {
    dplyr::select(out, "variant")
  } else {
    dplyr::select(out, "variant", dplyr::all_of(additional_columns))
  }
}
