library(variantstring)
library(rlang)
library(optparse)

options(dplyr.summarise.inform = FALSE)

opts <- list(
  make_option(
    c("-i", "--mlaf_input"),
    help = "TSV containing multi-locus allele frequency data with columns: group_id, variant, freq. Variant column is in STAVE format.",
    type = "character",
    default = NULL,
    callback = function(opt, flag_string, value, parser, ...) {
      if (!file.exists(value)) {
        stop(stringr::str_c(value, " does not exist"))
      }
      value
    }
  ),
  make_option(
    c("-o", "--output"),
    help = "Output file name. Default: %default",
    type = "character",
    default = "single_locus_allele_frequencies.tsv"
  )
)

args <- parse_args(OptionParser(option_list = opts))
required_args <- c("mlaf_input")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", paste(missing_args, collapse = ", ")))
}

load_mlaf <- function(path) {
  mlaf <- readr::read_tsv(
    path,
    col_types = readr::cols(
      group_id = readr::col_character(),
      variant = readr::col_character(),
      freq = readr::col_double()
    ),
    col_select = c("group_id", "variant", "freq")
  )
  return(mlaf)
}

convert_mlaf_to_slaf <- function(dat) {
  slaf <- dat |>
    dplyr::mutate(alleles = variantstring::variant_to_long(.data$variant)) |>
    tidyr::unnest("alleles") |>
    dplyr::group_by(.data$group_id, .data$gene, .data$pos, .data$aa) |>
    dplyr::summarize(freq = sum(.data$freq)) |>
    dplyr::ungroup() |>
    dplyr::arrange("group_id", "gene", "pos", "aa")
  return(slaf)
}

mlaf <- load_mlaf(args$mlaf_input)
slaf <- convert_mlaf_to_slaf(mlaf)
readr::write_tsv(slaf, args$output)
