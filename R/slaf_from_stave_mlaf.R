#' Load multi-locus allele frequency (MLAF) data
#'
#' @param path Path to a TSV with columns `group_id`, `variant`, and `freq`.
#' @return A tibble with those three columns.
#' @keywords internal
load_mlaf <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      group_id = readr::col_character(),
      variant = readr::col_character(),
      freq = readr::col_double()
    ),
    col_select = c("group_id", "variant", "freq"),
    show_col_types = FALSE
  )
}

#' Expand STAVE MLAF rows into single-locus allele frequencies
#'
#' Requires the optional **variantstring** package.
#'
#' @param dat MLAF tibble with `group_id`, `variant`, and `freq`.
#' @return Tibble with `group_id`, `gene_id`, `aa_position`, `aa`, and `freq`.
#' @keywords internal
convert_mlaf_to_slaf <- function(dat) {
  check_suggested_pkg(
    "variantstring",
    "expanding STAVE multi-locus variants via slaf_from_stave_mlaf()"
  )

  dat |>
    dplyr::mutate(alleles = variantstring::variant_to_long(.data$variant)) |>
    tidyr::unnest("alleles") |>
    dplyr::group_by(.data$group_id, .data$gene, .data$pos, .data$aa) |>
    dplyr::summarize(freq = sum(.data$freq), .groups = "drop") |>
    dplyr::rename(gene_id = "gene", aa_position = "pos") |>
    dplyr::arrange(
      .data$group_id, .data$gene_id, .data$aa_position, .data$aa
    )
}

#' Convert STAVE multi-locus allele frequencies to single-locus frequencies
#'
#' Expands STAVE `variant` strings with **variantstring**, aggregates
#' frequencies per amino acid allele, and emits STAVE-style single-locus
#' `variant` identifiers via [convert_single_locus_table_to_stave()].
#'
#' The **variantstring** package is an optional dependency (Suggests). It is
#' not installed automatically with PGEcore.
#'
#' @param mlaf Path to an MLAF TSV, or a data frame, with columns
#'   `group_id`, `variant`, and `freq`.
#' @param output Optional output TSV path. If `NULL`, results are returned
#'   without writing. CLI default is `single_locus_allele_frequencies.tsv`.
#'
#' @return A tibble with columns `variant` and `freq` (legacy STAVE conversion
#'   drops `group_id`, matching the original script).
#'
#' @export
slaf_from_stave_mlaf <- function(mlaf, output = NULL) {
  options(dplyr.summarise.inform = FALSE)
  check_suggested_pkg(
    "variantstring",
    "expanding STAVE multi-locus variants via slaf_from_stave_mlaf()"
  )

  if (is.character(mlaf) && length(mlaf) == 1L) {
    if (!file.exists(mlaf)) {
      stop(mlaf, " does not exist", call. = FALSE)
    }
    mlaf_tbl <- load_mlaf(mlaf)
  } else if (is.data.frame(mlaf)) {
    validate_required_columns(
      mlaf,
      c("group_id", "variant", "freq"),
      "mlaf"
    )
    mlaf_tbl <- tibble::as_tibble(mlaf)
  } else {
    stop("`mlaf` must be a file path or a data frame.", call. = FALSE)
  }

  slaf <- convert_mlaf_to_slaf(mlaf_tbl)
  # Match legacy script: only variant + freq (group_id not retained).
  slaf_output <- convert_single_locus_table_to_stave(
    slaf,
    additional_columns = "freq"
  )

  if (!is.null(output)) {
    readr::write_tsv(slaf_output, output)
  }
  slaf_output
}
