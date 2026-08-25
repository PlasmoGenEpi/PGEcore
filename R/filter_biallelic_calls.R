#' Read and validate amino acid call tables for biallelic filtering
#'
#' @param amino_acid_calls_fnp Path to a TSV with amino acid calls.
#' @return A tibble of amino acid calls.
#' @keywords internal
read_in_amino_acid_calls <- function(amino_acid_calls_fnp) {
  amino_acid_calls <- readr::read_tsv(amino_acid_calls_fnp, show_col_types = FALSE)
  validate_required_columns(
    amino_acid_calls,
    c("gene_id", "aa_position", "ref_aa", "aa"),
    amino_acid_calls_fnp
  )

  rules <- validate::validator(
    is.character(gene_id),
    is.numeric(aa_position),
    is.character(ref_aa),
    is.character(aa),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(ref_aa),
    !is.na(aa)
  )
  stop_on_validate_fails(amino_acid_calls, rules, "aa_calls")
  amino_acid_calls
}

#' Filter amino acid calls to biallelic loci
#'
#' Keeps loci (`gene_id`, `aa_position`, `ref_aa`) with at most two distinct
#' `aa` alleles. Optionally writes non-biallelic loci to a second file.
#'
#' @param aa_calls Path to a TSV, or a data frame, with columns
#'   `gene_id`, `aa_position`, `ref_aa`, and `aa`.
#' @param output Optional path for the biallelic output TSV.
#' @param nonbiallelic_output Optional path for loci with more than two alleles.
#' @param overwrite If `FALSE` (default), refuse to overwrite existing outputs.
#'
#' @return A list with tibbles `biallelic` and `nonbiallelic`. Each includes an
#'   `allele_calls` column with the distinct allele count per locus.
#'
#' @examples
#' path <- system.file("extdata", "example_aa_calls.tsv", package = "PGEcore")
#' filter_biallelic_calls(path)
#'
#' @export
filter_biallelic_calls <- function(aa_calls,
                                   output = NULL,
                                   nonbiallelic_output = NULL,
                                   overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)

  if (is.character(aa_calls) && length(aa_calls) == 1L) {
    if (!file.exists(aa_calls)) {
      stop(aa_calls, " does not exist", call. = FALSE)
    }
    aa_calls <- read_in_amino_acid_calls(aa_calls)
  } else if (is.data.frame(aa_calls)) {
    validate_required_columns(
      aa_calls,
      c("gene_id", "aa_position", "ref_aa", "aa"),
      "aa_calls"
    )
    aa_calls <- tibble::as_tibble(aa_calls)
  } else {
    stop("`aa_calls` must be a file path or a data frame.", call. = FALSE)
  }

  stop_if_output_exists(output, overwrite)
  stop_if_output_exists(nonbiallelic_output, overwrite)

  aa_calls <- aa_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position, .data$ref_aa) |>
    dplyr::mutate(allele_calls = dplyr::n_distinct(.data$aa)) |>
    dplyr::ungroup()

  biallelic <- dplyr::filter(aa_calls, .data$allele_calls <= 2)
  nonbiallelic <- dplyr::filter(aa_calls, .data$allele_calls > 2)

  if (!is.null(output)) {
    readr::write_tsv(biallelic, output)
  }
  if (!is.null(nonbiallelic_output)) {
    readr::write_tsv(nonbiallelic, nonbiallelic_output)
  }

  list(biallelic = biallelic, nonbiallelic = nonbiallelic)
}
