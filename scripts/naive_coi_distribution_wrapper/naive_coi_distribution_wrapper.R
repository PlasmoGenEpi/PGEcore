library(optparse)
library(rlang)

options(dplyr.summarise.inform = FALSE)

opts <- list(
  make_option(
    c("-i", "--coi_calls"),
    help = "TSV containing coi calls, with the columns: specimen_id, coi ",
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
    default = "coi_distribution.tsv"
  )
)

args <- parse_args(OptionParser(option_list = opts))
required_args <- c("coi_calls")
missing_args <- setdiff(required_args, names(args))
if (length(missing_args) > 0) {
  stop(stringr::str_c("Missing required arguments: ", missing_args))
}

load_coi_calls <- function(path) {
  coi_dat <- readr::read_tsv(
    path,
    col_types = readr::cols(
      specimen_id = readr::col_character(),
      coi = readr::col_integer()
    ),
    col_select = c("specimen_id", "coi")
  )
  return(coi_dat)
}

calculate_coi_distribution <- function(coi_calls) {
  coi_dist <- coi_calls |>
    dplyr::group_by(.data$coi) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::mutate(proportion = .data$n / sum(.data$n))
  return(coi_dist)
}

coi_calls <- load_coi_calls(args$coi_calls)
coi_dist <- calculate_coi_distribution(coi_calls)

readr::write_tsv(coi_dist, args$output)