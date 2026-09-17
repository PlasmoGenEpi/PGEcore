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

#' Collect validate::confront failures as a warning string (or NULL)
#'
#' Matches legacy scripts that recorded validation problems in `warns` and
#' continued rather than stopping.
#'
#' @param df Data frame that was confronted.
#' @param rules A [validate::validator()] object.
#' @param data_name Label used in the message.
#' @return Character warning message, or `NULL` if all rules pass.
#' @keywords internal
warn_on_validate_fails <- function(df, rules, data_name) {
  fails <- validate::confront(df, rules, raise = "all") |>
    validate::summary() |>
    dplyr::filter(.data$fails > 0)
  if (nrow(fails) > 0) {
    return(paste0(
      "Input ", data_name, " failed one or more validation checks: ",
      paste(fails$expression, collapse = "\n")
    ))
  }
  NULL
}

#' Stop when validate::confront reports failing rules
#'
#' @param df Data frame that was confronted.
#' @param rules A [validate::validator()] object.
#' @param data_name Label used in the error message.
#' @return Invisibly returns `TRUE` if all rules pass.
#' @keywords internal
stop_on_validate_fails <- function(df, rules, data_name) {
  msg <- warn_on_validate_fails(df, rules, data_name)
  if (!is.null(msg)) {
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}

#' Stop if an output path exists and overwrite is FALSE
#'
#' @param path Output file path (ignored when `NULL`).
#' @param overwrite Whether overwriting is allowed.
#' @return Invisibly returns `TRUE` if writing may proceed.
#' @keywords internal
stop_if_output_exists <- function(path, overwrite = FALSE) {
  if (!is.null(path) && nzchar(path) && file.exists(path) && !isTRUE(overwrite)) {
    stop(
      "file ", path, " already exists, use --overwrite to overwrite it",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Drop NULL entries from a parsed optparse argument list
#'
#' Used so [check_optparse_required_args()] can treat unset options as missing.
#'
#' @param arg Named list of parsed arguments.
#' @return `arg` with `NULL` values removed.
#' @keywords internal
drop_null_args <- function(arg) {
  arg[!vapply(arg, is.null, logical(1))]
}

#' Required columns for panel reference BED tables used when adding ref_seq
#' @keywords internal
ref_bed_required_cols <- function() {
  c("#chrom", "start", "end", "target_name", "length", "strand")
}

#' Read and validate a panel reference BED table
#'
#' @param path Path to a TSV with a header row.
#' @return A tibble of the BED table.
#' @keywords internal
read_ref_bed_table <- function(path) {
  ref_bed <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  validate_ref_bed_table(ref_bed, path)
  ref_bed
}

#' Validate panel reference BED columns and types
#'
#' @param ref_bed Data frame of panel locations.
#' @param data_name Label used in error messages.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_ref_bed_table <- function(ref_bed, data_name = "ref_bed") {
  validate_required_columns(ref_bed, ref_bed_required_cols(), data_name)
  rules <- validate::validator(
    is.character(`#chrom`),
    is.numeric(start),
    is.numeric(end),
    is.character(target_name),
    is.numeric(length),
    is.character(strand),
    !is.na(`#chrom`),
    !is.na(start),
    !is.na(end),
    !is.na(target_name),
    !is.na(length),
    !is.na(strand)
  )
  stop_on_validate_fails(ref_bed, rules, data_name)
  invisible(TRUE)
}

#' Load a FASTA with record names truncated at the first whitespace
#'
#' @param path Path to a FASTA file.
#' @param reason Passed to [check_suggested_pkg()] for **Biostrings**.
#' @return A `DNAStringSet` with shortened names.
#' @keywords internal
read_genome_dna_string_set <- function(path, reason) {
  check_suggested_pkg("Biostrings", reason)
  genome <- Biostrings::readDNAStringSet(path)
  names(genome) <- sub("\\s.*$", "", names(genome))
  genome
}

#' Stop if a character vector contains duplicated values
#'
#' @param values Values to check.
#' @param source Label for the input file or table.
#' @param field Field name used in the error message.
#' @return Invisibly returns `TRUE` if all values are unique.
#' @keywords internal
stop_on_duplicate_names <- function(values, source, field = "target_name") {
  dup_tbl <- tibble::tibble(value = values) |>
    dplyr::count(.data$value) |>
    dplyr::filter(.data$n > 1)
  if (nrow(dup_tbl) > 0) {
    stop(
      "found multi names for ", field, " in ", source,
      " found the following multiple times: ",
      paste(dup_tbl$value, collapse = ","),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Create an output directory, optionally replacing an existing one
#'
#' @param output_directory Directory path to create.
#' @param overwrite_dir If `TRUE`, delete `output_directory` first when it exists.
#' @return Invisibly returns `output_directory`.
#' @keywords internal
ensure_output_directory <- function(output_directory, overwrite_dir = FALSE) {
  if (dir.exists(output_directory) && isTRUE(overwrite_dir)) {
    unlink(output_directory, recursive = TRUE)
  } else if (dir.exists(output_directory)) {
    stop(
      output_directory,
      " already exist, use --overwrite_dir to overwrite",
      call. = FALSE
    )
  }
  dir.create(output_directory)
  invisible(output_directory)
}

#' Read a microhaplotype allele table (specimen_name, target_name, reads, seq)
#'
#' @param path Path to a TSV.
#' @return A tibble.
#' @keywords internal
read_mhap_allele_table <- function(path) {
  allele_table <- readr::read_tsv(
    path,
    col_types = readr::cols(specimen_name = readr::col_character()),
    show_col_types = FALSE
  )
  validate_required_columns(
    allele_table,
    c("specimen_name", "target_name", "reads", "seq"),
    path
  )
  allele_table
}

#' Validate a panel BED table that includes `ref_seq`
#'
#' @param ref_bed Data frame of panel locations with sequences.
#' @param data_name Label used in error messages.
#' @return Invisibly returns `TRUE` if valid.
#' @keywords internal
validate_ref_bed_with_seq_table <- function(ref_bed, data_name = "ref_bed") {
  validate_required_columns(
    ref_bed,
    c(ref_bed_required_cols(), "ref_seq"),
    data_name
  )
  invisible(TRUE)
}

#' Read a panel BED table that includes `ref_seq`
#'
#' @param path Path to a TSV with a header row.
#' @return A tibble of the BED table.
#' @keywords internal
read_ref_bed_with_seq_table <- function(path) {
  ref_bed <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  validate_ref_bed_with_seq_table(ref_bed, path)
  ref_bed
}
