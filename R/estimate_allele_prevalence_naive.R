#' Load amino acid calls for naive allele prevalence estimation
#'
#' @param path Path to a TSV of amino acid calls.
#' @return A tibble with columns `specimen_name`, `target_name`, and `variant`.
#' @keywords internal
prev_parse_aa_calls <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      gene_id = readr::col_character(),
      reads = readr::col_integer(),
      aa_position = readr::col_integer(),
      ref_aa = readr::col_character(),
      aa = readr::col_character()
    ),
    col_select = c(
      "specimen_name", "gene_id",
      "aa_position", "aa"
    )
  ) |>
    tidyr::unite("target_name", "gene_id", "aa_position", sep = ":") |>
    dplyr::rename(variant = "aa")
}

#' Load microhaplotype calls for naive allele prevalence estimation
#'
#' @param path Path to a TSV of microhaplotype genotypes.
#' @return A tibble with columns including `specimen_name`, `target_name`, and
#'   `variant`.
#' @keywords internal
prev_parse_mh_calls <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      target_name = readr::col_character(),
      seq = readr::col_character(),
      reads = readr::col_integer()
    )
  ) |>
    dplyr::rename(variant = "seq")
}

#' Calculate allele prevalence
#'
#' @param allele_table Allele table with `specimen_name`, `target_name`, and
#'   `variant`.
#' @return A tibble with `target_name`, `variant`, `prev`, `sample_count`, and
#'   `sample_total`.
#' @keywords internal
calculate_prevalence <- function(allele_table) {
  allele_table |>
    dplyr::group_by(.data$target_name) |>
    dplyr::mutate(sample_total = dplyr::n_distinct(.data$specimen_name)) |>
    dplyr::group_by(.data$target_name, .data$variant, .data$sample_total) |>
    dplyr::summarise(
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(prev = .data$sample_count / .data$sample_total) |>
    dplyr::relocate("prev", .before = "sample_total")
}

#' Format naive allele-prevalence output for AA or MH input
#'
#' @param prevalence Prevalence table with `target_name` and `variant`.
#' @param from_aa Logical; `TRUE` when input was amino acid calls.
#' @return Formatted tibble for writing.
#' @keywords internal
format_naive_prev_output <- function(prevalence, from_aa) {
  if (isTRUE(from_aa)) {
    prevalence |>
      tidyr::separate_wider_delim(
        "target_name",
        ":",
        names = c("gene_id", "aa_position")
      ) |>
      dplyr::rename(aa = "variant") |>
      convert_single_locus_table_to_stave(
        additional_columns = c("prev", "sample_count", "sample_total")
      )
  } else {
    dplyr::rename(prevalence, seq = "variant")
  }
}

#' Estimate allele prevalence naively from AA or microhaplotype calls
#'
#' Exactly one of `aa_calls` or `allele_table` must be provided. Prevalence is the
#' fraction of specimens carrying each allele at a locus.
#'
#' @param aa_calls Optional path to a TSV of amino acid calls with columns
#'   `specimen_name`, `gene_id`, `aa_position`, and `aa`.
#' @param allele_table Optional path to a TSV of microhaplotype genotypes with
#'   columns `specimen_name`, `target_name`, and `seq`.
#' @param output Optional output TSV path. If `NULL`, results are returned
#'   without writing a file.
#'
#' @return A tibble of allele prevalences. For amino acid input: `variant`,
#'   `prev`, `sample_count`, `sample_total`. For microhaplotype input:
#'   `target_name`, `seq`, `prev`, `sample_count`, `sample_total`.
#'
#' @examples
#' aa_path <- system.file(
#'   "extdata", "example_aa_calls.tsv",
#'   package = "PGEcore"
#' )
#' estimate_allele_prevalence_naive(aa_calls = aa_path)
#'
#' @export
estimate_allele_prevalence_naive <- function(aa_calls = NULL,
                                             allele_table = NULL,
                                             output = NULL) {
  options(dplyr.summarise.inform = FALSE)

  if (is.null(aa_calls) == is.null(allele_table)) {
    stop(
      "One and only one of the args --aa_calls and --allele_table should be ",
      "provided.",
      call. = FALSE
    )
  }

  from_aa <- !is.null(aa_calls)
  if (from_aa) {
    if (!file.exists(aa_calls)) {
      stop(aa_calls, " does not exist", call. = FALSE)
    }
    allele_table <- prev_parse_aa_calls(aa_calls)
  } else {
    if (!file.exists(allele_table)) {
      stop(allele_table, " does not exist", call. = FALSE)
    }
    allele_table <- prev_parse_mh_calls(allele_table)
  }

  prevalence <- calculate_prevalence(allele_table)
  prev_output <- format_naive_prev_output(prevalence, from_aa = from_aa)

  if (!is.null(output)) {
    readr::write_tsv(prev_output, output)
  }

  prev_output
}
