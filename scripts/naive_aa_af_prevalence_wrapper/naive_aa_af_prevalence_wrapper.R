library(optparse)
library(rlang)

opts <- list(
  make_option(
    "--aa_calls",
    help = stringr::str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ",
      "target_id, gene_id, aa_position, ref_codon, ref_aa, codon, aa"
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
    "--method",
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
    "--output",
    help = stringr::str_c(
      "Output file name. Default: %default"
    ),
    type = "character",
    default = "af_prevalence.tsv"
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
      .data$specimen_id, .data$gene_id, .data$aa_position, .data$ref_codon,
      .data$ref_aa
    ) |>
    # calculate the within-sample allele frequency
    dplyr::mutate(wsaf = .data$read_count / sum(.data$read_count)) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position, .data$ref_aa, .data$aa
    ) |>
    # sum af by gene and variant
    dplyr::summarise(af = sum(.data$wsaf)) |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    dplyr::mutate(af = .data$af / sum(.data$af)) |>
    dplyr::ungroup()
  return(af)
}

#' Calculate allele frequency by presence/absence
calculate_af_presence_absence <- function(aa_calls) {
  af <- aa_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    dplyr::mutate(total_obs = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position,
      .data$ref_aa, .data$aa, .data$total_obs
    ) |>
    dplyr::summarise(
      count = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(af = .data$count / .data$total_obs) |>
    dplyr::select(-"total_obs", -"count")

  return(af)
}

#' Calculate allele prevalence
calculate_prevalence <- function(aa_calls) {
  prev <- aa_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position, .data$ref_aa,) |>
    dplyr::mutate(n_samples = dplyr::n_distinct(.data$specimen_id)) |>
    dplyr::group_by(
      .data$gene_id, .data$aa_position, .data$ref_aa, .data$aa, .data$n_samples
    ) |>
    dplyr::summarise(
      count = dplyr::n()
    ) |>
    dplyr::mutate(prev = .data$count / .data$n_samples) |>
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
      ref_codon = readr::col_character(),
      ref_aa = readr::col_character(),
      codon = readr::col_character(),
      aa = readr::col_character()
    ),
    col_select = c(
      "specimen_id", "gene_id", "read_count",
      "aa_position", "ref_codon", "ref_aa", "codon", "aa"
    )
  )
  return(aa_dat)
}

aa_calls <- parse_aa_calls(args$aa_calls)

prevalence <- calculate_prevalence(aa_calls)
af <- switch(
  args$method,
  read_count_prop = calculate_af_read_count_prop(aa_calls),
  presence_absence = calculate_af_presence_absence(aa_calls)
)

out <- dplyr::left_join(
  prevalence, af,
  by = c("gene_id", "aa_position", "ref_aa", "aa")
)

readr::write_tsv(out, args$output)