#' Treat a missing CLI/R `NA` threshold as unset
#'
#' @param value Threshold or `NULL`.
#' @return `NULL` if `value` is `NA`; otherwise `value`.
#' @keywords internal
null_if_na_threshold <- function(value) {
  if (length(value) == 1L && is.na(value)) {
    return(NULL)
  }
  value
}

#' Validate `min_reads` / `min_wsaf` and treat `NA` as unset
#'
#' @param min_reads Minimum read count, or `NULL` / `NA` for no filter.
#' @param min_wsaf Minimum WSAF, or `NULL` / `NA` for no filter.
#' @return A list with `min_reads` and `min_wsaf` (`NULL` when unset).
#' @keywords internal
normalize_read_wsaf_thresholds <- function(min_reads, min_wsaf) {
  min_reads <- null_if_na_threshold(min_reads)
  min_wsaf <- null_if_na_threshold(min_wsaf)
  if (!is.null(min_reads) && (!is.numeric(min_reads) || min_reads < 0)) {
    stop("`min_reads` must be >= 0.", call. = FALSE)
  }
  if (!is.null(min_wsaf) && (!is.numeric(min_wsaf) || min_wsaf < 0 || min_wsaf > 1)) {
    stop("`min_wsaf` must be between 0 and 1 inclusive.", call. = FALSE)
  }
  list(min_reads = min_reads, min_wsaf = min_wsaf)
}

#' Drop rows below read-count and/or WSAF floors
#'
#' WSAF is `reads / sum(reads)` within `group_cols`. Thresholds are applied to
#' that original WSAF (not renormalized after dropping rows). Temporary
#' `total_reads` and `wsaf` columns are removed.
#'
#' @param df A tibble with a numeric `reads` column.
#' @param group_cols Character vector of columns that define a within-sample
#'   locus.
#' @param min_reads Minimum read count, or `NULL` for no filter.
#' @param min_wsaf Minimum WSAF, or `NULL` for no filter.
#' @param empty_message Error text if every row is dropped.
#' @return `df` with rows that failed a supplied threshold removed.
#' @keywords internal
filter_by_reads_and_wsaf <- function(df,
                                     group_cols,
                                     min_reads,
                                     min_wsaf,
                                     empty_message) {
  filtered <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::mutate(
      total_reads = sum(.data$reads, na.rm = TRUE),
      wsaf = dplyr::if_else(
        .data$total_reads > 0,
        .data$reads / .data$total_reads,
        NA_real_
      )
    ) |>
    dplyr::ungroup()

  if (!is.null(min_reads)) {
    filtered <- dplyr::filter(filtered, .data$reads >= min_reads)
  }
  if (!is.null(min_wsaf)) {
    filtered <- dplyr::filter(
      filtered,
      !is.na(.data$wsaf),
      .data$wsaf >= min_wsaf
    )
  }

  filtered <- dplyr::select(filtered, -c("total_reads", "wsaf"))

  if (nrow(filtered) == 0) {
    stop(empty_message, call. = FALSE)
  }
  filtered
}
