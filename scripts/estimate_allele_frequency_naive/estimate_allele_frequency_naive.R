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
      "target_id, gene_id, aa_position,, ref_aa, aa"
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
required_args <- c("aa_calls")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", missing_args))
}

#' Calculate allele frequency by within-sample allele frequency
calculate_af_read_count_prop <- function(aa_calls) {
  af <- aa_calls |>
    dplyr::group_by(
      .data$specimen_id, .data$gene_id, .data$aa_position, .data$ref_aa
    ) |>
    # calculate the within-sample allele frequency for each position
    dplyr::mutate(wsaf = .data$read_count / sum(.data$read_count)) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position, .data$ref_aa, .data$aa
    ) |>
    # sum af by gene and variant across all samples
    dplyr::summarise(freq = sum(.data$wsaf)) |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    # normalize by the total number of samples with data at each position
    dplyr::mutate(total = sum(.data$freq), freq = .data$freq / .data$total) |>
    dplyr::relocate("freq", .after = "total") |>
    dplyr::ungroup()
  return(af)
}

#' Calculate allele frequency by presence/absence
calculate_af_presence_absence <- function(aa_calls) {
  af <- aa_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    dplyr::mutate(total = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position,
      .data$ref_aa, .data$aa, .data$total
    ) |>
    dplyr::summarise(
      observed_alleles = dplyr::n(),
    ) |>
    dplyr::mutate(freq = .data$observed_alleles / .data$total) |>
    dplyr::select(-"observed_alleles")

  return(af)
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

out <- switch(
  args$method,
  read_count_prop = calculate_af_read_count_prop(aa_calls),
  presence_absence = calculate_af_presence_absence(aa_calls)
)
freq_output <- convert_single_locus_table_to_stave(out, "freq")
readr::write_tsv(freq_output, args$output)
