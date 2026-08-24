#' Create locus data from an allele table file
#'
#' Reads a TSV with columns `specimen_name`, `target_name`, and `seq`, renames
#' them to `sample_id`, `target_name`, and `allele`, and validates character
#' non-missing values with the **validate** package.
#'
#' @param input_path Path to the allele table TSV.
#' @return A data frame with columns `sample_id`, `target_name`, and `allele`.
#' @keywords internal
create_locus_data <- function(input_path) {
  print("Reading input data")
  input_data <- utils::read.csv(
    input_path,
    na.strings = "NA",
    sep = "\t",
    colClasses = c(specimen_name = "character")
  )
  locus_data <- input_data |>
    dplyr::select("specimen_name", "target_name", "seq") |>
    dplyr::rename(sample_id = "specimen_name", allele = "seq")

  print("Validating input format")
  rules <- validate::validator(
    is.character(sample_id),
    is.character(target_name),
    is.character(allele),
    !is.na(sample_id),
    !is.na(target_name),
    !is.na(allele)
  )

  print("Confronting input data with validation rules")
  fails <- validate::confront(locus_data, rules, raise = "all") |>
    validate::summary() |>
    dplyr::filter(.data$fails > 0)

  if (nrow(fails) > 0) {
    stop(
      "Analysis object failed one or more validation checks: ",
      paste(fails$expression, collapse = "\n"),
      call. = FALSE
    )
  }

  print("Returning Locus data")
  locus_data
}

#' Summarize allele metrics by target name
#'
#' @param locus_data A data frame with columns `target_name` and `allele`.
#' @return A tibble with `target_name`, `total_allele_count`,
#'   `unique_allele_count`, and `allele_singlets`.
#' @keywords internal
summarize_allele_table <- function(locus_data) {
  locus_data |>
    dplyr::group_by(.data$target_name) |>
    dplyr::summarize(
      total_allele_count = length(.data$allele),
      unique_allele_count = length(unique(.data$allele)),
      allele_singlets = sum(table(.data$allele) == 1)
    ) |>
    dplyr::ungroup()
}

#' Summarize alleles per locus from an allele table
#'
#' For each `target_name`, computes total allele count, unique allele count,
#' and the number of alleles that appear only once (singlets).
#'
#' @param allele_table Path to a TSV with columns `specimen_name`,
#'   `target_name`, and `seq`.
#' @param output Optional output TSV path. Defaults to
#'   `"allele_summary_by_target.tsv"` (matching the legacy CLI). If `NULL`,
#'   results are returned without writing a file.
#'
#' @return A tibble with columns `target_name`, `total_allele_count`,
#'   `unique_allele_count`, and `allele_singlets`.
#'
#' @examples
#' allele_path <- system.file(
#'   "extdata", "example_allele_table.tsv",
#'   package = "PGEcore"
#' )
#' allele_per_locus_summary(allele_path, output = NULL)
#'
#' @export
allele_per_locus_summary <- function(allele_table,
                                     output = "allele_summary_by_target.tsv") {
  if (!file.exists(allele_table)) {
    stop(allele_table, " does not exist", call. = FALSE)
  }

  locus_data <- create_locus_data(allele_table)
  allele_summary_by_target <- summarize_allele_table(locus_data)

  if (!is.null(output)) {
    readr::write_tsv(allele_summary_by_target, output)
  }

  allele_summary_by_target
}
