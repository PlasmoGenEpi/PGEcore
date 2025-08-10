# Get the script dir and make path to utils.R
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
utils_path <- file.path(script_dir, "..", "utils", "utils.R")

source(utils_path)
library(optparse)
library(rlang)

options(dplyr.summarise.inform = FALSE)

opts <- list(
  make_option(
    c("-i", "--aa_calls"),
    help = stringr::str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ",
      "target_id, gene_id, aa_position, ref_aa, aa"
    ),
    type = "character",
    default = NULL,
    callback = function(opt, flag_string, value, parser, ...) {
      if (!file.exists(value)) {
        stop(stringr::str_c(opt$aa_calls, " does not exist"))
      }
      value
    }
  ),
  make_option(
    c("-o", "--output"),
    help = stringr::str_c(
      "Output file name. Default: %default"
    ),
    type = "character",
    default = "prevalence.tsv"
  )
)

args <- parse_args(OptionParser(option_list = opts))
required_args <- c("aa_calls")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", missing_args))
}

#' Calculate allele prevalence
calculate_prevalence <- function(aa_calls) {
  prev <- aa_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position, .data$ref_aa) |>
    dplyr::mutate(total = dplyr::n_distinct(.data$specimen_id)) |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position, .data$ref_aa, .data$aa, .data$total
    ) |>
    dplyr::summarise(
      count = dplyr::n()
    ) |>
    dplyr::mutate(prev = .data$count / .data$total) |>
    dplyr::select(-"count")
  return(prev)
}

#' Load amino acid calls
parse_aa_calls <- function(path) {
  aa_dat <- readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_id = readr::col_character(),
      gene_id = readr::col_character(),
      read_count = readr::col_integer(),
      aa_position = readr::col_integer(),
      ref_aa = readr::col_character(),
      aa = readr::col_character()
    ),
    col_select = c(
      "specimen_id", "gene_id", "read_count",
      "aa_position", "ref_aa", "aa"
    )
  )
  return(aa_dat)
}


aa_calls <- parse_aa_calls(args$aa_calls)

prevalence <- calculate_prevalence(aa_calls)

prev_output <- convert_single_locus_table_to_stave(prevalence, "prev")

readr::write_tsv(prev_output, args$output)
