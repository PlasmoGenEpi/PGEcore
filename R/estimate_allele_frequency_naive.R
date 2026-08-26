#' Load amino acid calls for naive allele frequency estimation
#'
#' @param path Path to a TSV of amino acid calls.
#' @return A tibble with columns `specimen_name`, `target_name`, `reads`, and
#'   `variant`.
#' @keywords internal
af_parse_aa_calls <- function(path) {
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
      "specimen_name", "gene_id", "reads",
      "aa_position", "aa"
    )
  ) |>
    tidyr::unite("target_name", "gene_id", "aa_position", sep = ":") |>
    dplyr::rename(variant = "aa")
}

#' Load microhaplotype calls for naive allele frequency estimation
#'
#' @param path Path to a TSV of microhaplotype genotypes.
#' @return A tibble with columns including `specimen_name`, `target_name`,
#'   `reads`, and `variant`.
#' @keywords internal
af_parse_mh_calls <- function(path) {
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

#' Calculate allele frequency from within-sample allele proportions
#'
#' @param allele_table Allele table with `specimen_name`, `target_name`,
#'   `variant`, and `reads`.
#' @return A tibble with `target_name`, `variant`, and `freq`.
#' @keywords internal
calculate_af_read_count_prop <- function(allele_table) {
  allele_table |>
    dplyr::group_by(
      .data$specimen_name, .data$target_name
    ) |>
    dplyr::mutate(wsaf = .data$reads / sum(.data$reads)) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$target_name, .data$variant
    ) |>
    dplyr::summarise(freq = sum(.data$wsaf)) |>
    dplyr::group_by(.data$target_name) |>
    dplyr::mutate(total = sum(.data$freq), freq = .data$freq / .data$total) |>
    dplyr::select(-"total") |>
    dplyr::ungroup()
}

#' Calculate allele frequency from presence/absence
#'
#' @param allele_table Allele table with `target_name` and `variant`.
#' @return A tibble with `target_name`, `variant`, `allele_total`,
#'   `allele_count`, and `freq`.
#' @keywords internal
calculate_af_presence_absence <- function(allele_table) {
  allele_table |>
    dplyr::group_by(.data$target_name) |>
    dplyr::mutate(allele_total = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$target_name, .data$variant, .data$allele_total
    ) |>
    dplyr::summarise(
      allele_count = dplyr::n(),
    ) |>
    dplyr::mutate(freq = .data$allele_count / .data$allele_total)
}

#' Format naive allele-frequency output for AA or MH input
#'
#' @param out Frequency table with `target_name` and `variant`.
#' @param from_aa Logical; `TRUE` when input was amino acid calls.
#' @return Formatted tibble for writing.
#' @keywords internal
format_naive_af_output <- function(out, from_aa) {
  if (isTRUE(from_aa)) {
    out |>
      tidyr::separate_wider_delim(
        "target_name",
        ":",
        names = c("gene_id", "aa_position")
      ) |>
      dplyr::rename(aa = "variant") |>
      tidyr::unite("variant", "gene_id", "aa_position", "aa", sep = ":")
  } else {
    dplyr::rename(out, seq = "variant")
  }
}

#' Estimate allele frequency naively from AA or microhaplotype calls
#'
#' Exactly one of `aa_calls` or `allele_table` must be provided. Frequency is
#' estimated either from within-sample read-count proportions
#' (`read_count_prop`) or from presence/absence (`presence_absence`).
#'
#' ## Inputs
#'
#' - **`aa_calls`** (optional): AA calls (`specimen_name`, `gene_id`,
#'   `aa_position`, `aa`, `reads`). See
#'   `vignette("input-formats", package = "PGEcore")`.
#' - **`allele_table`** (optional): Allele table (`specimen_name`,
#'   `target_name`, `seq`, `reads`). See the same vignette.
#'
#' ## Outputs
#'
#' - **`output`** (optional): Allele-frequency TSV. For AA input, `variant` is
#'   a STAVE-style `gene_id:aa_position:aa` string plus `freq` (and count
#'   columns for `presence_absence`). For microhaplotype input: `target_name`,
#'   `seq`, `freq` (plus counts for `presence_absence`). If `NULL`, results are
#'   returned without writing a file.
#'
#' ## Running
#'
#' ```r
#' estimate_allele_frequency_naive(
#'   aa_calls = "aa_calls.tsv",
#'   output = "allele_frequency.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/estimate_allele_frequency_naive \
#'   --aa_calls aa_calls.tsv \
#'   --output allele_frequency.tsv
#' ```
#'
#' @param aa_calls Optional path to an AA calls TSV. See *Inputs*.
#' @param allele_table Optional path to an allele table TSV. See *Inputs*.
#' @param method Estimation method: `"presence_absence"` (default) or
#'   `"read_count_prop"`.
#' @param output Optional output TSV path. Default for the CLI is
#'   `allele_frequency.tsv`.
#'
#' @return A tibble of estimated allele frequencies. For amino acid input the
#'   `variant` column is a STAVE-style `gene_id:aa_position:aa` string. For
#'   microhaplotype input columns include `target_name`, `seq`, and `freq`
#'   (plus count columns for `presence_absence`).
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @examples
#' aa_path <- system.file(
#'   "extdata", "example_aa_calls.tsv",
#'   package = "PGEcore"
#' )
#' estimate_allele_frequency_naive(aa_calls = aa_path)
#'
#' @export
estimate_allele_frequency_naive <- function(aa_calls = NULL,
                                            allele_table = NULL,
                                            method = "presence_absence",
                                            output = NULL) {
  options(dplyr.summarise.inform = FALSE)

  if (is.null(aa_calls) == is.null(allele_table)) {
    stop(
      "One and only one of the args --aa_calls and --allele_table should be ",
      "provided.",
      call. = FALSE
    )
  }
  if (!method %in% c("read_count_prop", "presence_absence")) {
    stop(method, " is not a valid method", call. = FALSE)
  }

  from_aa <- !is.null(aa_calls)
  if (from_aa) {
    if (!file.exists(aa_calls)) {
      stop(aa_calls, " does not exist", call. = FALSE)
    }
    allele_table <- af_parse_aa_calls(aa_calls)
  } else {
    if (!file.exists(allele_table)) {
      stop(allele_table, " does not exist", call. = FALSE)
    }
    allele_table <- af_parse_mh_calls(allele_table)
  }

  out <- switch(
    method,
    read_count_prop = calculate_af_read_count_prop(allele_table),
    presence_absence = calculate_af_presence_absence(allele_table)
  )

  freq_output <- format_naive_af_output(out, from_aa = from_aa)

  if (!is.null(output)) {
    readr::write_tsv(freq_output, output)
  }

  freq_output
}
