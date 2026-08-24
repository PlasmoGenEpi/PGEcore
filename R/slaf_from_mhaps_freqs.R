#' Read and validate microhaplotype allele frequency tables
#'
#' @param mhaps_slaf_fnp Path to a TSV with columns `target_name`, `seq`,
#'   `freq`, and `sample_total`.
#' @return A tibble of microhaplotype allele frequencies.
#' @keywords internal
process_input_mhaps_slaf <- function(mhaps_slaf_fnp) {
  input <- readr::read_tsv(mhaps_slaf_fnp, show_col_types = FALSE)
  validate_required_columns(
    input,
    c("target_name", "seq", "freq", "sample_total"),
    mhaps_slaf_fnp
  )
  rules <- validate::validator(
    is.character(target_name),
    is.character(seq),
    is.numeric(freq),
    is.numeric(sample_total),
    !is.na(target_name),
    !is.na(seq),
    !is.na(freq),
    !is.na(sample_total)
  )
  warns <- warn_on_validate_fails(input, rules, mhaps_slaf_fnp)
  if (!is.null(warns)) {
    warning(warns, call. = FALSE)
  }
  input
}

#' Read and validate translated loci-of-interest for microhaplotypes
#'
#' @param loci_of_interest_per_microhaps_fnp Path to a TSV with columns
#'   `target_name`, `gene_id`, `aa_position`, `seq`, and `aa`.
#' @return A tibble of translated loci.
#' @keywords internal
process_input_loci_of_interest_per_microhaps <- function(loci_of_interest_per_microhaps_fnp) {
  input <- readr::read_tsv(
    loci_of_interest_per_microhaps_fnp,
    show_col_types = FALSE
  )
  validate_required_columns(
    input,
    c("target_name", "gene_id", "aa_position", "seq", "aa"),
    loci_of_interest_per_microhaps_fnp
  )
  rules <- validate::validator(
    is.character(target_name),
    is.character(gene_id),
    is.numeric(aa_position),
    is.character(seq),
    is.character(aa),
    !is.na(target_name),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(seq),
    !is.na(aa)
  )
  warns <- warn_on_validate_fails(input, rules, loci_of_interest_per_microhaps_fnp)
  if (!is.null(warns)) {
    warning(warns, call. = FALSE)
  }
  input
}

#' Calculate single-locus allele frequencies from microhaplotype frequencies
#'
#' Joins microhaplotype allele frequencies to translated amino acid calls,
#' aggregates frequencies per amino acid allele, and renormalises so frequencies
#' sum to one. Collapsed output averages evenly across overlapping targets
#' (taking `max(sample_total)`).
#'
#' @param mhaps_slaf Path or data frame with columns `target_name`, `seq`,
#'   `freq`, and `sample_total`.
#' @param loci_of_interest_per_microhaps Path or data frame with columns
#'   `target_name`, `gene_id`, `aa_position`, `seq`, and `aa`.
#' @param slaf_output Optional path for collapsed SLAF TSV (`variant`, `freq`,
#'   `sample_total`).
#' @param per_target_slaf_output Optional path for per-target SLAF TSV.
#'
#' @return A list with tibbles `slaf` (collapsed) and `per_target_slaf`.
#'
#' @export
slaf_from_mhaps_freqs <- function(mhaps_slaf,
                                  loci_of_interest_per_microhaps,
                                  slaf_output = NULL,
                                  per_target_slaf_output = NULL) {
  options(dplyr.summarise.inform = FALSE)

  if (is.character(mhaps_slaf) && length(mhaps_slaf) == 1L) {
    if (!file.exists(mhaps_slaf)) {
      stop(mhaps_slaf, " does not exist", call. = FALSE)
    }
    mhaps_tbl <- process_input_mhaps_slaf(mhaps_slaf)
  } else if (is.data.frame(mhaps_slaf)) {
    validate_required_columns(
      mhaps_slaf,
      c("target_name", "seq", "freq", "sample_total"),
      "mhaps_slaf"
    )
    mhaps_tbl <- tibble::as_tibble(mhaps_slaf)
  } else {
    stop("`mhaps_slaf` must be a file path or a data frame.", call. = FALSE)
  }

  if (
    is.character(loci_of_interest_per_microhaps) &&
      length(loci_of_interest_per_microhaps) == 1L
  ) {
    if (!file.exists(loci_of_interest_per_microhaps)) {
      stop(loci_of_interest_per_microhaps, " does not exist", call. = FALSE)
    }
    loci_tbl <- process_input_loci_of_interest_per_microhaps(
      loci_of_interest_per_microhaps
    )
  } else if (is.data.frame(loci_of_interest_per_microhaps)) {
    validate_required_columns(
      loci_of_interest_per_microhaps,
      c("target_name", "gene_id", "aa_position", "seq", "aa"),
      "loci_of_interest_per_microhaps"
    )
    loci_tbl <- tibble::as_tibble(loci_of_interest_per_microhaps)
  } else {
    stop(
      "`loci_of_interest_per_microhaps` must be a file path or a data frame.",
      call. = FALSE
    )
  }

  combined_tables <- mhaps_tbl |>
    dplyr::inner_join(loci_tbl, by = c("target_name", "seq"))

  per_target_slaf <- combined_tables |>
    dplyr::group_by(
      .data$target_name, .data$gene_id, .data$aa_position,
      .data$sample_total, .data$aa
    ) |>
    dplyr::summarise(freq = sum(.data$freq), .groups = "drop") |>
    dplyr::group_by(
      .data$target_name, .data$gene_id, .data$aa_position, .data$sample_total
    ) |>
    dplyr::mutate(
      total_freq = sum(.data$freq),
      freq = .data$freq / .data$total_freq
    ) |>
    dplyr::select(-"total_freq") |>
    tidyr::unite("variant", "gene_id", "aa_position", "aa", sep = ":") |>
    dplyr::ungroup()

  # Even weighting across overlapping targets; take max sample_total.
  slaf <- combined_tables |>
    dplyr::group_by(.data$gene_id, .data$aa_position, .data$aa) |>
    dplyr::summarise(
      freq = sum(.data$freq),
      sample_total = max(.data$sample_total),
      .groups = "drop"
    ) |>
    dplyr::group_by(.data$gene_id, .data$aa_position, .data$sample_total) |>
    dplyr::mutate(
      total_freq = sum(.data$freq),
      freq = .data$freq / .data$total_freq
    ) |>
    dplyr::select(-"total_freq") |>
    tidyr::unite("variant", "gene_id", "aa_position", "aa", sep = ":") |>
    dplyr::ungroup()

  if (!is.null(slaf_output)) {
    readr::write_tsv(slaf, slaf_output)
  }
  if (!is.null(per_target_slaf_output)) {
    readr::write_tsv(per_target_slaf, per_target_slaf_output)
  }

  list(slaf = slaf, per_target_slaf = per_target_slaf)
}
