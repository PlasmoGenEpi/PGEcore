
library(stringr)
library(optparse)
library(dplyr)
library(validate)

# ----------------------------------------------------------------------
# Parse arguments

opts <- list(
  make_option(
    "--input_path",
    type = "character",
    help = str_c(
      "Path to a TSV file containing allele calls, with the columns: specimen_id, ", 
      "target_id, read_count, seq"
    )
  ), 
  make_option(
    "--output_path",
    type = "character",
    help = "Path to write a TSV file containing results, with the columns: specimin_id, coi"
  ), 
  make_option(
    "--method",
    type = "character",
    default = "integer_method",
    help = str_c(
      "Which method is used to infer COI: 'integer_method' or 'quantile_method'. ",
      "Both methods start by calculating the number of observed alleles at each ",
      "locus, and sorting these in decreasing order. The integer method then takes ",
      "the nth value as the COI estimate. A higher n will tend to lead to lower COI ",
      "estimates, but is more robust to errors. The quantile method works in the same ",
      "way, but by specifying a percentage quantile rather than a raw integer. This has ",
      "the advantage of naturally scaling with the number of loci, while the integer method ",
      "does not."
      ),
    metavar = "method",
  ),
  make_option(
    "--integer_threshold",
    type = "integer",
    default = 1,
    help = str_c("The index in the ordered sequence that is used as the COI estimate ",
                 "(used in integer_method only). A positive integer of 1 or above."),
    metavar = "integer"
  ),
  make_option(
    "--quantile_threshold",
    type = "double",
    default = 0.05,
    help = str_c("The quantile used to define the COI estimate from the ordered sequence ",
                 "(used in quantile_method only). A numeric value between 0 and 1"),
    metavar = "number"
  )
)

# ----------------------------------------------------------------------
#' Estimate COI using naive methods
#'
#' Takes the path of TSV file of allele calls. For every sample, counts the
#' number of alleles observed at each locus and sorts them in decreasing order.
#' Two methods are implemented to estimate COI from these values:
#' `method = integer_method` and `method = quantile_method`. The integer method takes the
#' nth value as the COI estimate, meaning higher values will tend to lead to lower COI
#' estimates but are more robust to errors. In the quantile method, a
#' percentage quantile is used rather than a raw integer. This has the advantage
#' of naturally accounting for the number of loci.
#'
#' @param input_path Path to a TSV file with allele calls. It should 
#'   have columns for specimen_id, target_id, read_count, seq.
#' @param method Which method to use. One of `integer_method` or `quantile_method`.
#' @param integer_threshold Which sorted value to use as the COI estimate. Only
#'   used if `method = integer_method`.
#' @param quantile_threshold The quantile (between 0 and 1) used to make the COI
#'   estimate. Only used if `method = quantile_method`.
#' 
#' @return data.frame of COI estimates for each sample.

run_estimate_coi_naive <- function(input_path,
                                   method = "integer_method",
                                   integer_threshold = 1,
                                   quantile_threshold = 0.05) {
  
  # Check input arguments
  stopifnot(is.character(input_path))
  stopifnot(method %in% c("integer_method", "quantile_method"))
  stopifnot(is.numeric(integer_threshold))
  stopifnot(integer_threshold == as.integer(integer_threshold))
  stopifnot(integer_threshold > 0)
  stopifnot(is.numeric(quantile_threshold))
  stopifnot((quantile_threshold >= 0) & (quantile_threshold <= 1))
  
  # Read allele calls from file
  df_alleles <- read.table(file = input_path, header = TRUE)
  
  # Validate input format
  rules <- validate::validator(is.character(specimen_id),
                               is.character(target_id),
                               is.numeric(read_count) & read_count == as.integer(read_count) & read_count > 0,
                               is.character(seq),
                               ! is.na(specimen_id),
                               ! is.na(target_id),
                               ! is.na(read_count),
                               ! is.na(seq)
                               )
  df_fails <- validate::confront(df_alleles, rules, raise = "all") |>
    validate::summary()
  
  if (any(df_fails$fails)) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(df_fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  
  # Count number of alleles at each locus
  df_n_alleles <- df_alleles |>
    dplyr::group_by(specimen_id, target_id) |>
    dplyr::summarise(n_alleles = n_distinct(seq), .groups = "drop")
  
  # Get number of loci per sample, and work out the integer corresponding to the
  # desired quantile (may or may not be used)
  df_loci <- df_alleles |>
    dplyr::group_by(specimen_id) |>
    dplyr::summarise(loci = n()) |>
    dplyr::mutate(n_limit = floor((loci - 1)*quantile_threshold) + 1)
  
  if (any(integer_threshold > df_loci$loci)) {
    stop("integer_threshold exceeds number of loci for one or more samples")
  }
  
  # Estimate COI by one or other method
  if (method == "integer_method") {
    
    df_coi <- df_n_alleles |>
      dplyr::group_by(specimen_id) |>
      dplyr::summarise(coi = sort(n_alleles, decreasing = TRUE)[integer_threshold])
    
  } else if (method == "quantile_method") {
    
    df_coi <- df_n_alleles |>
      dplyr::group_by(specimen_id) |>
      dplyr::arrange(desc(n_alleles)) |>
      dplyr::mutate(row_number = row_number()) |>
      dplyr::left_join(df_loci, by = join_by(specimen_id)) |>
      dplyr::filter(row_number == n_limit) |>
      dplyr::rename(coi = n_alleles) |>
      dplyr::select(specimen_id, coi)
  }
  
  return(df_coi)
}

# ----------------------------------------------------------------------
# Check that an object matches the expected format output from run_estimate_coi_naive()
check_coi_format <- function(df_coi) {
  stopifnot(is.data.frame(df_coi))
  stopifnot(ncol(df_coi) == 2)
  stopifnot(all(colnames(df_coi) == c("specimen_id", "coi")))
  stopifnot(!any(is.na(df_coi$specimen_id)))
  stopifnot(all(is.character(df_coi$specimen_id)))
  stopifnot(!any(is.na(df_coi$coi)))
  stopifnot(all(is.numeric(df_coi$coi)))
  stopifnot(all(df_coi$coi == as.integer(df_coi$coi)))
  stopifnot(all(df_coi$coi > 0))
  invisible(TRUE)
}

# ----------------------------------------------------------------------
#' Write COI estimates to file
#'
#' Takes a data.frame of COI estimates and writes to TSV.
#'
#' @param df_coi a data.frame of estimated COI values. Must contain the columns
#'   `specimen_id` and `coi`.
#' @param output_path file path to write output.
#' 
#' @return data.frame of COI estimates for each sample.

write_coi_naive <- function(df_coi,
                            output_path) {
  
  # Check inputs
  check_coi_format(df_coi)
  stopifnot(is.character(output_path))
  
  # Write to file as TSV
  write.table(x = df_coi,
              file = output_path,
              row.names = FALSE)
}

#------------------------------------------------------
# RUN MODULE

# parse input arguments
args <- parse_args(OptionParser(option_list = opts))

# estimate COI
df_coi <- run_estimate_coi_naive(input_path = args$input_path,
                                 method = args$method,
                                 integer_threshold = args$integer_threshold,
                                 quantile_threshold = args$quantile_threshold)

# Write to file
write_coi_naive(df_coi = df_coi,
                output_path = args$output_path)


