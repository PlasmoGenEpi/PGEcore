#' Validate amino-acid call columns and types for filtering
#'
#' @param aa_calls A data frame of amino-acid calls.
#' @param data_name Label used in error messages.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_aa_calls_for_filter <- function(aa_calls, data_name) {
  validate_required_columns(
    aa_calls,
    c("specimen_name", "gene_id", "aa_position", "aa", "reads"),
    data_name
  )
  rules <- validate::validator(
    is.character(specimen_name) || is.factor(specimen_name),
    is.character(gene_id) || is.factor(gene_id),
    is.numeric(aa_position),
    is.character(aa) || is.factor(aa),
    is.numeric(reads),
    !is.na(specimen_name),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(aa),
    !is.na(reads),
    reads >= 0
  )
  stop_on_validate_fails(aa_calls, rules, data_name)
  invisible(TRUE)
}

#' Filter amino-acid calls by read count and/or within-sample allele frequency
#'
#' Drops amino-acid calls below a minimum read count (`min_reads`) and/or a
#' minimum within-sample allele frequency (`min_wsaf`). WSAF is the call's
#' `reads` divided by the total `reads` for that `specimen_name`, `gene_id`,
#' and `aa_position`. When both thresholds are set, a call must pass both.
#' When neither is set, the table is returned unchanged. Extra columns are
#' kept.
#'
#' ## Inputs
#'
#' - **`aa_calls`**: AA calls (`specimen_name`, `gene_id`, `aa_position`, `aa`,
#'   `reads`), as a file path or data frame. See
#'   `vignette("input-formats", package = "PGEcore")`.
#'
#' ## Outputs
#'
#' - **`output`** (optional): Filtered AA calls TSV (same columns as the
#'   input).
#'
#' ## Running
#'
#' ```r
#' filter_aa_calls(
#'   aa_calls = "aa_calls.tsv",
#'   output = "filtered_aa_calls.tsv",
#'   min_reads = 20,
#'   min_wsaf = 0.05
#' )
#' ```
#'
#' ```bash
#' Rscript exec/filter_aa_calls \
#'   --aa_calls aa_calls.tsv \
#'   --output filtered_aa_calls.tsv \
#'   --min_reads 20 \
#'   --min_wsaf 0.05
#' ```
#'
#' @param aa_calls Path to an AA calls TSV, or a data frame with the same
#'   columns. See *Inputs*.
#' @param output Optional path for the filtered AA calls TSV.
#' @param min_reads Minimum read count to keep a call. `NULL` (default): no
#'   read-count filter.
#' @param min_wsaf Minimum within-sample allele frequency (0-1) to keep a
#'   call. `NULL` (default): no WSAF filter.
#' @param overwrite If `FALSE` (default), refuse to overwrite an existing
#'   `output`.
#'
#' @return A tibble of amino-acid calls that passed the supplied thresholds.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @examples
#' path <- system.file("extdata", "example_aa_calls.tsv", package = "PGEcore")
#' filter_aa_calls(path, min_reads = 3)
#'
#' @export
filter_aa_calls <- function(aa_calls,
                            output = NULL,
                            min_reads = NULL,
                            min_wsaf = NULL,
                            overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  thresholds <- normalize_read_wsaf_thresholds(min_reads, min_wsaf)
  aa_calls <- as_input_tibble(
    aa_calls,
    function(path) readr::read_tsv(path, show_col_types = FALSE),
    "aa_calls"
  )
  validate_aa_calls_for_filter(aa_calls, "aa_calls")
  stop_if_output_exists(output, overwrite)

  filtered <- filter_by_reads_and_wsaf(
    aa_calls,
    group_cols = c("specimen_name", "gene_id", "aa_position"),
    min_reads = thresholds$min_reads,
    min_wsaf = thresholds$min_wsaf,
    empty_message = "All amino-acid calls were removed by the read-count / WSAF filters."
  )

  if (!is.null(output)) {
    readr::write_tsv(filtered, output)
  }
  filtered
}
