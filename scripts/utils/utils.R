suppressPackageStartupMessages(library(dplyr))

#' Convert prev df to STAVE Format for variant column
#'
#' This function processes a data frame by generating a `variant` column, 
#' which combines `gene_id`, `aa_position`, and `aa` into a single identifier.
#' It then returns only the `variant` and `prev` columns.
#'
#' @param df A data frame containing at least the columns:
#'   - `gene_id` (character): Gene identifier
#'   - `aa_position` (numeric or integer): Amino acid position
#'   - `aa` (character): Amino acid at the given position
#'   - `prev` (numeric): Prevalence of the variant
#'
#' @return A data frame with two columns:
#'   - `variant` (character): Concatenated string of `gene_id:aa_position:aa`
#'   - `prev` (numeric): The prevalence value corresponding to the variant
#'
#' @examples
#' df <- data.frame(
#'   gene_id = c("PF3D7_0417200.1", "PF3D7_0417200.1"),
#'   aa_position = c(51, 59),
#'   aa = c("I", "R"),
#'   prev = c(0.5, 1)
#' )
#' convert_to_stave(df)
convert_prev_table_to_stave <- function(df) {

  df %>%
    ungroup %>%
    mutate(variant = paste(gene_id, aa_position, aa, sep = ":")) %>%
    select(variant, prev)

}