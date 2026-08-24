#' Read and validate a locus table for per-locus popgen summaries
#'
#' @param input_path Path to a TSV allele table.
#' @param specimen_name_col Specimen ID column name.
#' @param target_name_col Target/locus column name.
#' @param target_value_col Allele/sequence column name.
#' @return A tibble with `specimen_name`, `target_name`, and `target_value`.
#' @keywords internal
read_popgen_locus_data <- function(input_path,
                                   specimen_name_col = "specimen_name",
                                   target_name_col = "target_name",
                                   target_value_col = "seq") {
  input_data <- readr::read_tsv(
    input_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      reads = readr::col_double()
    ),
    progress = FALSE,
    show_col_types = FALSE
  )
  needed <- c(specimen_name_col, target_name_col, target_value_col)
  validate_required_columns(input_data, needed, input_path)
  locus_data <- input_data |>
    dplyr::select(dplyr::all_of(needed)) |>
    dplyr::rename_with(
      ~ c("specimen_name", "target_name", "target_value"),
      .cols = dplyr::all_of(needed)
    )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(target_value),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(target_value)
  )
  stop_on_validate_fails(locus_data, rules, "locus_data")
  locus_data
}

#' External MSA binary expected for a given `msa::msa()` method
#'
#' @param msa_method One of `"Muscle"`, `"ClustalW"`, or `"ClustalOmega"`.
#' @return Executable name for [check_external_tool()].
#' @keywords internal
msa_method_binary <- function(msa_method) {
  switch(
    msa_method,
    Muscle = "muscle",
    ClustalW = "clustalw",
    ClustalOmega = "clustalo",
    stop(
      "--msa_method must be 'ClustalW', 'ClustalOmega', or 'Muscle', not ",
      msa_method,
      call. = FALSE
    )
  )
}

#' Population-genetic statistics for one locus
#'
#' Aligns unique allele sequences with **msa**, restores duplicate sequences,
#' then computes nucleotide diversity, segregating sites, and Tajima's D with
#' **pegas**. Identical sequences skip alignment and return zeros.
#'
#' Requires the **msa** alignment binary on `PATH` (`muscle`, `clustalw`, or
#' `clustalo` depending on `msa_method`). The **msa** R package does not bundle
#' those tools; install them separately (for example via Conda).
#'
#' @param allele_data Character vector of allele sequences.
#' @param msa_method Passed to [msa::msa()]: `"Muscle"`, `"ClustalW"`, or
#'   `"ClustalOmega"`.
#' @return A named list of statistics.
#' @keywords internal
calculate_popgen_stats <- function(allele_data, msa_method = "Muscle") {
  unique_seqs <- !duplicated(allele_data)
  unique_ids <- which(unique_seqs)
  allele_data_unique <- allele_data[unique_ids]
  names(allele_data_unique) <- seq_along(allele_data_unique)
  orig2unique_mapping <- match(allele_data, allele_data_unique)

  if (length(allele_data_unique) == 1) {
    return(list(
      Nucleotide_Diversity = 0,
      Segregating_Sites = 0,
      Tajima_D = 0
    ))
  }

  check_external_tool(msa_method_binary(msa_method))

  aligned_unique <- msa::msa(
    allele_data_unique,
    method = msa_method,
    type = "dna",
    order = "input"
  ) |>
    msa::msaConvert(type = "ape::DNAbin")
  if (!identical(names(allele_data_unique), labels(aligned_unique))) {
    stop("Order does not match between alignment input and output", call. = FALSE)
  }
  aligned_all <- aligned_unique[orig2unique_mapping, ]

  nucleotide_diversity <- pegas::nuc.div(aligned_all)
  segregating_sites <- length(ape::seg.sites(aligned_all))
  tajima_test <- pegas::tajima.test(aligned_all)

  list(
    Nucleotide_Diversity = nucleotide_diversity,
    Segregating_Sites = segregating_sites,
    Tajima_D = tajima_test$D,
    Tajima_D_pval_normal = tajima_test$Pval.normal,
    Tajima_D_pval_beta = tajima_test$Pval.beta
  )
}

#' Population-genetic statistics grouped by `target_name`
#'
#' @param locus_data Tibble with `target_name` and `target_value`.
#' @param msa_method Alignment method for `calculate_popgen_stats()`.
#' @return A tibble of per-locus statistics.
#' @keywords internal
calculate_stats_by_target_name <- function(locus_data, msa_method = "Muscle") {
  locus_data |>
    dplyr::group_by(.data$target_name) |>
    dplyr::summarise(
      stats = list(calculate_popgen_stats(.data$target_value, msa_method))
    ) |>
    tidyr::unnest_wider("stats")
}

#' Per-locus nucleotide diversity, segregating sites, and Tajima's D
#'
#' Groups allele sequences by locus and computes population-genetic summaries.
#' Requires **ape**, **msa**, and **pegas** (Suggests). **msa** calls an
#' external aligner; install the binary for `msa_method` on `PATH`:
#'
#' * `"Muscle"` — `muscle`
#' * `"ClustalW"` — `clustalw`
#' * `"ClustalOmega"` — `clustalo`
#'
#' @param allele_table Path to a TSV of alleles.
#' @param specimen_name_col Specimen ID column name.
#' @param target_name_col Locus column name.
#' @param target_value_col Allele/sequence column name.
#' @param out Optional output TSV path. Defaults to
#'   `"per_locus_popgen_summary.tsv"`. If `NULL`, results are returned without
#'   writing.
#' @param msa_method Alignment method: `"Muscle"` (default), `"ClustalW"`, or
#'   `"ClustalOmega"`.
#'
#' @return A tibble of per-locus statistics with lower-case column names.
#'
#' @export
per_locus_popgen_summary <- function(allele_table,
                                     specimen_name_col = "specimen_name",
                                     target_name_col = "target_name",
                                     target_value_col = "seq",
                                     out = "per_locus_popgen_summary.tsv",
                                     msa_method = "Muscle") {
  options(dplyr.summarise.inform = FALSE)
  if (!(msa_method %in% c("ClustalW", "ClustalOmega", "Muscle"))) {
    stop(
      "--msa_method must be 'ClustalW', 'ClustalOmega', or 'Muscle', not ",
      msa_method,
      call. = FALSE
    )
  }
  check_suggested_pkg("ape", "per_locus_popgen_summary()")
  check_suggested_pkg("msa", "per_locus_popgen_summary()")
  check_suggested_pkg("pegas", "per_locus_popgen_summary()")

  locus_data <- read_popgen_locus_data(
    allele_table,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col
  )
  res <- calculate_stats_by_target_name(locus_data, msa_method = msa_method)
  colnames(res) <- tolower(colnames(res))
  if (!is.null(out)) {
    readr::write_tsv(res, out)
  }
  res
}
