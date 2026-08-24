#' Check that an object matches the expected COI estimate format
#'
#' @param df_coi A data frame of COI estimates.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
check_coi_format <- function(df_coi) {
  stopifnot(is.data.frame(df_coi))
  stopifnot(ncol(df_coi) == 2)
  stopifnot(all(colnames(df_coi) == c("specimen_name", "coi")))
  stopifnot(!any(is.na(df_coi$specimen_name)))
  stopifnot(all(is.character(df_coi$specimen_name)))
  stopifnot(!any(is.na(df_coi$coi)))
  stopifnot(all(is.numeric(df_coi$coi)))
  stopifnot(all(df_coi$coi == as.integer(df_coi$coi)))
  stopifnot(all(df_coi$coi > 0))
  invisible(TRUE)
}

#' Validate allele-call columns used for naive COI estimation
#'
#' @param df_alleles Allele-call data frame.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_allele_calls_for_coi_naive <- function(df_alleles) {
  validate_required_columns(
    df_alleles,
    c("specimen_name", "target_name", "reads", "seq"),
    "allele calls"
  )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.numeric(reads) & reads == as.integer(reads) & reads > 0,
    is.character(seq),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(reads),
    !is.na(seq)
  )
  df_fails <- validate::confront(df_alleles, rules, raise = "all") |>
    validate::summary()

  if (any(df_fails$fails)) {
    stop(
      "Input input_data failed one or more validation checks: ",
      paste(df_fails$expression, collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Estimate COI from an allele-call data frame using naive methods
#'
#' @param df_alleles Data frame with columns `specimen_name`, `target_name`,
#'   `reads`, and `seq`.
#' @param method One of `integer_method` or `quantile_method`.
#' @param integer_threshold Index into the decreasing allele-count sequence
#'   (integer method only).
#' @param quantile_threshold Quantile in `[0, 1]` (quantile method only).
#' @return A tibble with columns `specimen_name` and `coi`.
#' @keywords internal
estimate_coi_naive_from_alleles <- function(df_alleles,
                                            method = "integer_method",
                                            integer_threshold = 1,
                                            quantile_threshold = 0.05) {
  stopifnot(method %in% c("integer_method", "quantile_method"))
  stopifnot(is.numeric(integer_threshold))
  stopifnot(integer_threshold == as.integer(integer_threshold))
  stopifnot(integer_threshold > 0)
  stopifnot(is.numeric(quantile_threshold))
  stopifnot((quantile_threshold >= 0) & (quantile_threshold <= 1))

  validate_allele_calls_for_coi_naive(df_alleles)

  df_n_alleles <- df_alleles |>
    dplyr::group_by(.data$specimen_name, .data$target_name) |>
    dplyr::summarise(
      n_alleles = dplyr::n_distinct(.data$seq),
      .groups = "drop"
    )

  # Preserve legacy behaviour: loci is allele-row count per specimen, not
  # unique target count.
  df_loci <- df_alleles |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::summarise(loci = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      n_limit = floor((.data$loci - 1) * quantile_threshold) + 1
    )

  if (any(integer_threshold > df_loci$loci)) {
    stop(
      "integer_threshold exceeds number of loci for one or more samples",
      call. = FALSE
    )
  }

  if (method == "integer_method") {
    df_coi <- df_n_alleles |>
      dplyr::group_by(.data$specimen_name) |>
      dplyr::summarise(
        coi = sort(.data$n_alleles, decreasing = TRUE)[integer_threshold],
        .groups = "drop"
      )
  } else {
    df_coi <- df_n_alleles |>
      dplyr::group_by(.data$specimen_name) |>
      dplyr::arrange(dplyr::desc(.data$n_alleles)) |>
      dplyr::mutate(row_number = dplyr::row_number()) |>
      dplyr::left_join(df_loci, by = "specimen_name") |>
      dplyr::filter(.data$row_number == .data$n_limit) |>
      dplyr::rename(coi = "n_alleles") |>
      dplyr::select("specimen_name", "coi")
  }

  df_coi
}

#' Estimate COI using naive allele-count methods
#'
#' For every specimen, counts distinct alleles at each locus and sorts those
#' counts in decreasing order. With `method = "integer_method"`, the
#' `integer_threshold`-th value is the COI estimate. With
#' `method = "quantile_method"`, the value at `quantile_threshold` is used
#' instead (scaling naturally with the number of observed allele rows).
#'
#' @param input_path Path to a TSV of allele calls with columns
#'   `specimen_name`, `target_name`, `reads`, and `seq`.
#' @param output_path Optional path to write a TSV with columns
#'   `specimen_name` and `coi`. If `NULL`, results are returned without writing.
#' @param method One of `"integer_method"` or `"quantile_method"`.
#'   Default: `"integer_method"`.
#' @param integer_threshold Positive integer index into the ordered allele
#'   counts (integer method only). Default: `1`.
#' @param quantile_threshold Quantile in `[0, 1]` (quantile method only).
#'   Values near zero yield higher COI estimates. Default: `0.05`.
#'
#' @return A tibble with columns `specimen_name` and `coi`.
#'
#' @examples
#' allele_path <- system.file(
#'   "extdata", "example_allele_table.tsv",
#'   package = "PGEcore"
#' )
#' estimate_coi_naive(allele_path, method = "integer_method")
#'
#' @export
estimate_coi_naive <- function(input_path,
                               output_path = NULL,
                               method = "integer_method",
                               integer_threshold = 1,
                               quantile_threshold = 0.05) {
  stopifnot(is.character(input_path), length(input_path) == 1L)
  if (!file.exists(input_path)) {
    stop(input_path, " does not exist", call. = FALSE)
  }

  df_alleles <- readr::read_tsv(
    input_path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      .default = readr::col_character(),
      reads = readr::col_integer()
    ),
    progress = FALSE
  )

  df_coi <- estimate_coi_naive_from_alleles(
    df_alleles = df_alleles,
    method = method,
    integer_threshold = integer_threshold,
    quantile_threshold = quantile_threshold
  )

  if (!is.null(output_path)) {
    check_coi_format(df_coi)
    stopifnot(is.character(output_path), length(output_path) == 1L)
    readr::write_tsv(df_coi, output_path)
  }

  df_coi
}
