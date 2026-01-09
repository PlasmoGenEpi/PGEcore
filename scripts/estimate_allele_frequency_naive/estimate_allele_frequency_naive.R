# Get the script dir and make path to utils.R
script_dir <- dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))
utils_path <- file.path(script_dir, "..", "utils", "utils.R")

source(utils_path)
library(optparse)
library(rlang)
library(tidyr)

options(dplyr.summarise.inform = FALSE)

opts <- list(
  make_option(
    c("--aa_calls"),
    help = stringr::str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ",
      "target_id, gene_id, aa_position, ref_aa, aa, read_count"
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
    c("--mh_calls"),
    help = stringr::str_c(
      "TSV containing microhaplotype genotypes, with the columns: ", 
      "specimen_id, target_id, seq, read_count"
    ),
    type = "character",
    default = NULL,
    callback = function(opt, flag_string, value, parser, ...) {
      if (!file.exists(value)) {
        stop(stringr::str_c(opt$mh_calls, " does not exist"))
      }
      value
    }
  ),
  make_option(
    c("-m", "--method"),
    help = stringr::str_c(
      "Method to use for estimating allele frequency. Options are: ",
      "read_count_prop, presence_absence. Default: %default"
    ),
    type = "character",
    default = "read_count_prop",
    callback = function(opt, flag_string, value, parser, ...) {
      if (!value %in% c("read_count_prop", "presence_absence")) {
        stop(stringr::str_c(value, " is not a valid method"))
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
    default = "allele_frequency.tsv"
  )
)

args <- parse_args(OptionParser(option_list = opts))
input_data_args <- c("aa_calls", "mh_calls")
absent_data_args <- setdiff(input_data_args, names(args))
if (length(absent_data_args) != 1) {
  stop(
    "One and only one of the args --aa_calls and --mh_calls should be ", 
    "provided."
  )
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
        "aa_position", "aa"
      )
    ) |>
    unite(target_id, gene_id, aa_position, sep = ":") |>
    rename(variant = aa)
  return(aa_dat)
}

#' Load microhaplotype calls
parse_mh_calls <- function(path) {
  mh_dat <- readr::read_tsv(
      path,
      col_types = readr::cols(
        specimen_id = readr::col_character(),
        target_id = readr::col_character(),
        seq = readr::col_character(),
        read_count = readr::col_integer()
      )
    ) |>
    dplyr::rename(variant = seq)
  return(mh_dat)
}

#' Calculate allele frequency by within-sample allele frequency
calculate_af_read_count_prop <- function(allele_table) {
  af <- allele_table |>
    dplyr::group_by(
      .data$specimen_id, .data$target_id
    ) |>
    # calculate the within-sample allele frequency for each position
    dplyr::mutate(wsaf = .data$read_count / sum(.data$read_count)) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$target_id, .data$variant
    ) |>
    # sum af by gene and variant across all samples
    dplyr::summarise(freq = sum(.data$wsaf)) |>
    dplyr::group_by(.data$target_id) |>
    # normalize by the total number of samples with data at each position
    dplyr::mutate(total = sum(.data$freq), freq = .data$freq / .data$total) |>
    dplyr::select(-total) |>
    dplyr::ungroup()
  return(af)
}

#' Calculate allele frequency by presence/absence
calculate_af_presence_absence <- function(allele_table) {
  af <- allele_table |>
    dplyr::group_by(.data$target_id) |>
    dplyr::mutate(total = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$target_id, .data$variant, .data$total
    ) |>
    dplyr::summarise(
      observed_alleles = dplyr::n(),
    ) |>
    dplyr::mutate(freq = .data$observed_alleles / .data$total) |>
    dplyr::select(-"observed_alleles")

  return(af)
}

# Read input
if (! is.null(args$aa_calls)) {
  allele_table <- parse_aa_calls(args$aa_calls)
} else {
  allele_table <- parse_mh_calls(args$mh_calls)
}

# Estimate allele frequency
out <- switch(
  args$method,
  read_count_prop = calculate_af_read_count_prop(allele_table),
  presence_absence = calculate_af_presence_absence(allele_table)
)

# Format and write output
if (! is.null(args$aa_calls)) {
  freq_output <- out |>
    tidyr::separate_wider_delim(
      target_id, 
      ":", 
      names = c("gene_id", "aa_position")
    ) |>
    dplyr::rename(aa = variant) |>
    convert_single_locus_table_to_stave("freq")
} else {
  freq_output <- out |>
    dplyr::rename(seq = variant)
}
readr::write_tsv(freq_output, args$output)
