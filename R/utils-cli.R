#' Stop if required optparse arguments are missing
#'
#' @param arg Named list of parsed arguments (as from [optparse::parse_args()]).
#' @param required_args Character vector of required argument names (without `--`).
#' @return Invisibly returns `TRUE` if all required arguments are present.
#' @keywords internal
check_optparse_required_args <- function(arg, required_args) {
  missing <- setdiff(required_args, names(arg))
  if (length(missing) > 0) {
    missing_flags <- paste0("--", missing)
    stop(
      "Missing the following arguments: ",
      paste(missing_flags, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Return names of required columns that are missing from a data frame
#'
#' @param df A data frame.
#' @param required_cols Character vector of required column names.
#' @return Character vector of missing column names (possibly empty).
#' @keywords internal
return_missing_columns <- function(df, required_cols) {
  setdiff(required_cols, colnames(df))
}

#' Validate that a data frame has required columns and is non-empty
#'
#' @param data Data frame to validate.
#' @param required_cols Character vector of required column names.
#' @param data_name Label used in error messages.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_required_columns <- function(data, required_cols, data_name) {
  missing_cols <- return_missing_columns(data, required_cols)
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", data_name, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(data) == 0) {
    stop(data_name, " is empty", call. = FALSE)
  }
  invisible(TRUE)
}

#' Decompose shared and unique values between two vectors
#'
#' @param vector_a First vector.
#' @param vector_b Second vector.
#' @return A list with `only_in_vector_a`, `only_in_vector_b`, `shared`, and `all`.
#' @keywords internal
set_decompose <- function(vector_a, vector_b) {
  list(
    only_in_vector_a = setdiff(vector_a, vector_b),
    only_in_vector_b = setdiff(vector_b, vector_a),
    shared = intersect(vector_a, vector_b),
    all = union(vector_a, vector_b)
  )
}
