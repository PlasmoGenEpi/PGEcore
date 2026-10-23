#' Validate allele table columns and types for filtering
#'
#' @param allele_table A data frame of allele calls.
#' @param data_name Label used in error messages.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_allele_table_for_filter <- function(allele_table, data_name) {
  validate_required_columns(
    allele_table,
    c("specimen_name", "target_name", "seq", "reads"),
    data_name
  )
  rules <- validate::validator(
    is.character(specimen_name) || is.factor(specimen_name),
    is.character(target_name) || is.factor(target_name),
    is.character(seq) || is.factor(seq),
    is.numeric(reads),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(seq),
    !is.na(reads),
    reads >= 0
  )
  stop_on_validate_fails(allele_table, rules, data_name)
  invisible(TRUE)
}

#' Filter an allele table by read count and/or within-sample allele frequency
#'
#' Drops alleles below a minimum read count (`min_reads`) and/or a minimum
#' within-sample allele frequency (`min_wsaf`). WSAF is the allele's `reads`
#' divided by the total `reads` for that `specimen_name` and `target_name`.
#' When both thresholds are set, an allele must pass both. When neither is set,
#' the table is returned unchanged. Extra columns are kept.
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele table (`specimen_name`, `target_name`, `seq`,
#'   `reads`), as a file path or data frame. See
#'   `vignette("input-formats", package = "PGEcore")`.
#'
#' ## Outputs
#'
#' - **`output`** (optional): Filtered allele table TSV (same columns as the
#'   input).
#'
#' ## Running
#'
#' ```r
#' filter_allele_table(
#'   allele_table = "allele_table.tsv",
#'   output = "filtered_allele_table.tsv",
#'   min_reads = 20,
#'   min_wsaf = 0.05
#' )
#' ```
#'
#' ```bash
#' Rscript exec/filter_allele_table \
#'   --allele_table allele_table.tsv \
#'   --output filtered_allele_table.tsv \
#'   --min_reads 20 \
#'   --min_wsaf 0.05
#' ```
#'
#' @param allele_table Path to an allele table TSV, or a data frame with the
#'   same columns. See *Inputs*.
#' @param output Optional path for the filtered allele table TSV.
#' @param min_reads Minimum read count to keep an allele. `NULL` (default): no
#'   read-count filter.
#' @param min_wsaf Minimum within-sample allele frequency (0-1) to keep an
#'   allele. `NULL` (default): no WSAF filter.
#' @param overwrite If `FALSE` (default), refuse to overwrite an existing
#'   `output`.
#'
#' @return A tibble of alleles that passed the supplied thresholds.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @examples
#' path <- system.file("extdata", "example_allele_table.tsv", package = "PGEcore")
#' filter_allele_table(path, min_reads = 50)
#'
#' @export
filter_allele_table <- function(allele_table,
                                output = NULL,
                                min_reads = NULL,
                                min_wsaf = NULL,
                                overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  thresholds <- normalize_read_wsaf_thresholds(min_reads, min_wsaf)
  allele_table <- as_input_tibble(
    allele_table,
    read_mhap_allele_table,
    "allele_table"
  )
  validate_allele_table_for_filter(allele_table, "allele_table")
  stop_if_output_exists(output, overwrite)

  filtered <- filter_by_reads_and_wsaf(
    allele_table,
    group_cols = c("specimen_name", "target_name"),
    min_reads = thresholds$min_reads,
    min_wsaf = thresholds$min_wsaf,
    empty_message = "All alleles were removed by the read-count / WSAF filters."
  )

  if (!is.null(output)) {
    readr::write_tsv(filtered, output)
  }
  filtered
}
