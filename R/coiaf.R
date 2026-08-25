#' Estimate complexity of infection (COI) using coiaf
#'
#' Processes SNP read-count data and optionally population-level minor allele
#' frequencies (PLMAF), then estimates COI with both the frequency and variant
#' methods from the **coiaf** package.
#'
#' The **coiaf** package is an optional dependency (Suggests). It is not
#' installed automatically with PGEcore.
#'
#' @param snp_calls A data frame with columns `specimen_name`, `snp_name`,
#'   `reads`, and `seq_base`.
#' @param plmaf Optional data frame with columns `snp_name`, `seq_base`, and
#'   `plmaf`. If `NULL`, PLMAF is calculated from `snp_calls`.
#' @param seq_error Sequencing error rate (default: `0.01`).
#' @param max_coi Maximum COI to consider (default: `25`).
#'
#' @return A data frame with columns `specimen_name`, `coi_freq`, and
#'   `coi_variant`.
#'
#' @export
run_coiaf <- function(snp_calls, plmaf = NULL, seq_error = 0.01, max_coi = 25) {
  check_suggested_pkg("coiaf", "COI estimation via run_coiaf()")

  validate_required_columns(
    snp_calls,
    c("specimen_name", "snp_name", "reads", "seq_base"),
    "SNP data"
  )
  if (!is.null(plmaf)) {
    validate_required_columns(
      plmaf,
      c("snp_name", "seq_base", "plmaf"),
      "PLMAF data"
    )
  }
  if (seq_error < 0 || seq_error > 1) {
    stop("seq_error must be between 0 and 1", call. = FALSE)
  }
  if (max_coi < 1) {
    stop("max_coi must be at least 1", call. = FALSE)
  }

  message("Processing SNP data...")

  complete_snp_data <- snp_calls |>
    tidyr::complete(
      specimen_name,
      tidyr::nesting(snp_name, seq_base),
      fill = list(reads = 0)
    ) |>
    dplyr::select("specimen_name", "snp_name", "seq_base", "reads") |>
    dplyr::group_by(.data$specimen_name, .data$snp_name) |>
    dplyr::mutate(n_snps = dplyr::n())

  if (any(complete_snp_data$n_snps > 2)) {
    warning(
      "Some targets have more than 2 SNPs. These will be excluded from the analysis.",
      call. = FALSE
    )
  }
  if (any(complete_snp_data$n_snps < 2)) {
    warning(
      "Some targets have fewer than 2 SNPs. These will be excluded from the analysis.",
      call. = FALSE
    )
  }

  complete_snp_data <- complete_snp_data |>
    dplyr::filter(.data$n_snps == 2) |>
    dplyr::select(-"n_snps")

  processed <- complete_snp_data |>
    dplyr::group_by(.data$specimen_name, .data$snp_name) |>
    dplyr::mutate(coverage = sum(.data$reads)) |>
    dplyr::ungroup() |>
    dplyr::mutate(wsmaf = .data$reads / .data$coverage) |>
    dplyr::select("specimen_name", "snp_name", "seq_base", "wsmaf", "reads")

  if (is.null(plmaf)) {
    message("Calculating population-level minor allele frequencies...")
    plmaf <- complete_snp_data |>
      dplyr::group_by(.data$snp_name, .data$seq_base) |>
      dplyr::summarize(reads = sum(.data$reads), .groups = "drop") |>
      dplyr::group_by(.data$snp_name) |>
      dplyr::mutate(plmaf = .data$reads / sum(.data$reads)) |>
      dplyr::arrange(.data$plmaf, .by_group = TRUE) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select("snp_name", "seq_base", "plmaf")
  }

  plmaf_snp_names <- unique(plmaf$snp_name)
  processed_snp_names <- unique(processed$snp_name)

  missing_in_plmaf <- setdiff(processed_snp_names, plmaf_snp_names)
  if (length(missing_in_plmaf) > 0) {
    warning(
      "Some SNPs in the observed data are not present in the PLMAF data, ",
      "these will be excluded from the analysis: ",
      paste(missing_in_plmaf, collapse = ", "),
      call. = FALSE
    )
  }

  missing_in_obs <- setdiff(plmaf_snp_names, processed_snp_names)
  if (length(missing_in_obs) > 0) {
    warning(
      "Some SNPs in the PLMAF data are not present in the observed data, ",
      "these will be excluded from the analysis: ",
      paste(missing_in_obs, collapse = ", "),
      call. = FALSE
    )
  }

  filtered_processed <- processed |>
    dplyr::left_join(plmaf, by = c("snp_name", "seq_base")) |>
    dplyr::filter(!is.na(.data$plmaf))

  message(
    "Estimating COI for ",
    length(unique(filtered_processed$specimen_name)),
    " specimens..."
  )

  filtered_processed |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::summarize(
      coi_freq = coiaf::optimize_coi(
        tibble::tibble(wsmaf, plmaf, coverage = reads),
        data_type = "real",
        coi_method = "frequency",
        seq_error = seq_error,
        max_coi = max_coi
      ),
      coi_variant = coiaf::optimize_coi(
        tibble::tibble(wsmaf, plmaf, coverage = reads),
        data_type = "real",
        coi_method = "variant",
        seq_error = seq_error,
        max_coi = max_coi
      ),
      .groups = "drop"
    )
}

#' Run COIAF from input and output file paths
#'
#' File-oriented entry point used by the `coiaf_wrapper` CLI.
#'
#' @param snp_calls Path to SNP data TSV.
#' @param output Path for output TSV.
#' @param plmaf Optional path to PLMAF TSV.
#' @param seq_error Sequencing error rate (default: `0.01`).
#' @param max_coi Maximum COI to consider (default: `25`).
#'
#' @return The result tibble (also written to `output`).
#' @export
coiaf_wrapper <- function(snp_calls,
                          output,
                          plmaf = NULL,
                          seq_error = 0.01,
                          max_coi = 25) {
  if (!file.exists(snp_calls)) {
    stop("SNP data file not found: ", snp_calls, call. = FALSE)
  }
  if (!is.null(plmaf) && !file.exists(plmaf)) {
    stop("PLMAF file not found: ", plmaf, call. = FALSE)
  }

  output_dir <- dirname(output)
  if (!identical(output_dir, ".") && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  snp_tbl <- readr::read_tsv(
    snp_calls,
    show_col_types = FALSE,
    col_types = readr::cols(specimen_name = readr::col_character())
  )

  plmaf_tbl <- NULL
  if (!is.null(plmaf)) {
    plmaf_tbl <- readr::read_tsv(plmaf, show_col_types = FALSE)
  }

  results <- run_coiaf(
    snp_calls = snp_tbl,
    plmaf = plmaf_tbl,
    seq_error = seq_error,
    max_coi = max_coi
  )

  readr::write_tsv(results, output)
  results
}
