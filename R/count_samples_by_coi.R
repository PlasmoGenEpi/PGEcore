#' Load COI calls from a TSV file
#'
#' Reads `specimen_name` and `coi`, and rounds `coi` to the nearest integer.
#'
#' @param path Path to a TSV with columns `specimen_name` and `coi`.
#' @return A tibble with columns `specimen_name` and integer-rounded `coi`.
#' @keywords internal
load_coi_table <- function(path) {
  coi_dat <- readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      coi = readr::col_double()
    ),
    col_select = c("specimen_name", "coi")
  )
  coi_dat$coi <- round(coi_dat$coi)
  coi_dat
}

#' Calculate the distribution of COI values across specimens
#'
#' @param coi_table A data frame with a numeric `coi` column.
#' @return A tibble with columns `coi`, `n`, and `proportion` for each integer
#'   COI from 1 to `max(coi)`.
#' @keywords internal
calculate_coi_distribution <- function(coi_table) {
  ret <- tibble::tibble(coi = seq_len(max(coi_table$coi)))
  coi_table |>
    dplyr::group_by(.data$coi) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
    dplyr::arrange(.data$coi) |>
    dplyr::right_join(ret, by = "coi") |>
    tidyr::replace_na(list(proportion = 0, n = 0)) |>
    dplyr::arrange(.data$coi)
}

#' Count specimens by complexity of infection (COI)
#'
#' Reads per-specimen COI values, rounds them to integers, and returns the
#' count and proportion of specimens at each COI level.
#'
#' ## Inputs
#'
#' - **`coi_table`**: COI table (`specimen_name`, `coi`), as a file path or
#'   data frame. See `vignette("input-formats", package = "PGEcore")`.
#'
#' ## Outputs
#'
#' - **`output`** (optional): TSV with columns `coi`, `n`, and `proportion`.
#'   If `NULL`, results are returned without writing a file.
#'
#' ## Running
#'
#' ```r
#' count_samples_by_coi(
#'   coi_table = "coi_table.tsv",
#'   output = "coi_distribution.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/count_samples_by_coi \
#'   --coi_table coi_table.tsv \
#'   --output coi_distribution.tsv
#' ```
#'
#' @param coi_table Path to a COI table TSV, or a data frame with the same
#'   columns. See *Inputs*.
#' @param output Optional output TSV path. Default for the CLI is
#'   `coi_distribution.tsv`.
#'
#' @return A tibble with columns `coi`, `n`, and `proportion`.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @examples
#' coi_path <- system.file("extdata", "example_coi_table.tsv", package = "PGEcore")
#' count_samples_by_coi(coi_path)
#'
#' @export
count_samples_by_coi <- function(coi_table, output = NULL) {
  options(dplyr.summarise.inform = FALSE)

  if (is.character(coi_table) && length(coi_table) == 1L) {
    if (!file.exists(coi_table)) {
      stop(coi_table, " does not exist", call. = FALSE)
    }
    coi_dat <- load_coi_table(coi_table)
  } else if (is.data.frame(coi_table)) {
    validate_required_columns(
      coi_table,
      c("specimen_name", "coi"),
      "COI table"
    )
    coi_dat <- tibble::as_tibble(coi_table)
    coi_dat$coi <- round(as.numeric(coi_dat$coi))
  } else {
    stop("`coi_table` must be a file path or a data frame.", call. = FALSE)
  }

  coi_dist <- calculate_coi_distribution(coi_dat)

  if (!is.null(output)) {
    readr::write_tsv(coi_dist, output)
  }

  coi_dist
}
