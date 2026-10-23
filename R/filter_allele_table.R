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

  if (length(min_reads) == 1L && is.na(min_reads)) {
    min_reads <- NULL
  }
  if (length(min_wsaf) == 1L && is.na(min_wsaf)) {
    min_wsaf <- NULL
  }

  if (!is.null(min_reads) && (!is.numeric(min_reads) || min_reads < 0)) {
    stop("`min_reads` must be >= 0.", call. = FALSE)
  }
  if (!is.null(min_wsaf) && (!is.numeric(min_wsaf) || min_wsaf < 0 || min_wsaf > 1)) {
    stop("`min_wsaf` must be between 0 and 1 inclusive.", call. = FALSE)
  }

  if (is.character(allele_table) && length(allele_table) == 1L) {
    if (!file.exists(allele_table)) {
      stop(allele_table, " does not exist", call. = FALSE)
    }
    allele_table <- read_mhap_allele_table(allele_table)
  } else if (is.data.frame(allele_table)) {
    allele_table <- tibble::as_tibble(allele_table)
  } else {
    stop("`allele_table` must be a file path or a data frame.", call. = FALSE)
  }
  validate_allele_table_for_filter(allele_table, "allele_table")

  stop_if_output_exists(output, overwrite)

  filtered <- allele_table |>
    dplyr::group_by(.data$specimen_name, .data$target_name) |>
    dplyr::mutate(
      total_reads = sum(.data$reads, na.rm = TRUE),
      wsaf = dplyr::if_else(
        .data$total_reads > 0,
        .data$reads / .data$total_reads,
        NA_real_
      )
    ) |>
    dplyr::ungroup()

  if (!is.null(min_reads)) {
    filtered <- dplyr::filter(filtered, .data$reads >= min_reads)
  }
  if (!is.null(min_wsaf)) {
    filtered <- dplyr::filter(
      filtered,
      !is.na(.data$wsaf),
      .data$wsaf >= min_wsaf
    )
  }

  filtered <- dplyr::select(filtered, -c("total_reads", "wsaf"))

  if (nrow(filtered) == 0) {
    stop(
      "All alleles were removed by the read-count / WSAF filters.",
      call. = FALSE
    )
  }

  if (!is.null(output)) {
    readr::write_tsv(filtered, output)
  }

  filtered
}
