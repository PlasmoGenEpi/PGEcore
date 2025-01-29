library(dplyr)
library(optparse)
library(readr)
library(stringr)
library(validate)

# install a specific version of variantstring package (this is in development so
# may not always be backwards-compatible)
if (!requireNamespace("variantstring", quietly = TRUE) ||
    packageVersion("variantstring") != "1.8.0") {
  remotes::install_github("mrc-ide/variantstring@1.8.0")
}
library(variantstring) 

# ----------------------------------------------------------------------
# Parse arguments
opts <- list(
  make_option(
    "--input_path", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing amino acid calls, with the columns: specimen_id, ", 
      "gene, pos, read_count, aa"
    )
  ), 
  make_option(
    "--output_path", 
    type = "character",
    help = str_c(
      "Path to write an output TSV file containing prevalence and frequency estimates ",
      "for all variants in the data"
    )
  )
)

# ----------------------------------------------------------------------
#' Convert long form amino acid table to variant strings
#'
#' @description
#' Takes the path of a TSV file of amino acid calls. Reads these in and converts
#' them to \href{https://github.com/mrc-ide/variantstring}{variant strings},
#' which are then returned.
#'
#' @param input_path Path to a TSV file with amino acid calls. It should have
#'   columns for specimen_id, gene, pos, read_count, aa.
#' 
#' @import dplyr
#' 
#' @return vector of variant strings corresponding to each specimen.

aa_table_to_variant <- function(input_path) {
  
  # Check input arguments
  stopifnot(is.character(input_path))
  
  # read in amino acid calls and validate columns
  df_aa <- read.table(input_path, header = TRUE)
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(gene_id), 
    is.integer(aa_position), 
    is.integer(read_count), 
    is.character(aa), 
    ! is.na(specimen_id), 
    ! is.na(gene_id), 
    ! is.na(aa_position), 
    ! is.na(read_count), 
    ! is.na(aa)
  )
  fails <- validate::confront(df_aa, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  # tidy up columns
  df_aa <- df_aa |>
    select(specimen_id, gene_id, aa_position, read_count, aa) |>
    rename(gene = gene_id,
           pos = aa_position)
  
  # get numer of amino acids at each locus
  df_n_aa <- df_aa |>
    group_by(specimen_id, gene, pos) |>
    summarise(n_aa = n(),
              .groups = "drop")
  
  # get final list of long format data.frames
  long_list <- df_aa |>
    left_join(df_n_aa, by = c("specimen_id", "gene", "pos")) |>
    mutate(het = (n_aa > 1),
           phased = FALSE) |>
    select(gene, pos, n_aa, het, phased, aa, read_count) |>
    split(f = df_aa$specimen)
  names(long_list) <- NULL
  
  # convert to variant strings
  vs <- long_to_variant(long_list)
  
  return(vs)
}

# ----------------------------------------------------------------------
#' Get prevalence of a set of variants
#'
#' @description
#' Takes a vector of target variants for which we want to calculate prevalence.
#' Compares each against a dataset, also in variant string format. Calculates
#' both the prevalence and the allele frequency of each target.
#'
#' @param target_variants a vector of target variant strings.
#' @param comparison_variants a vector of variant strings representing a full dataset.
#' 
#' @return data.frame of prevalence and frequency.

calculate_variant_prevalence <- function(target_variants,
                                         comparison_variants) {
  
  # Check input arguments
  check_variant_string(target_variants)
  check_variant_string(comparison_variants)
  
  # match each target against the comparison data
  l <- list()
  for (i in seq_along(target_variants)) {
    
    # make data.frame recording if this target is a match against all comparison
    # variants
    df_match <- compare_variant_string(target_string = target_variants[i],
                                       comparison_strings = comparison_variants)
    df_match$match_pos <- compare_position_string(target_string = position_from_variant_string(target_variants[i]),
                                                  comparison_strings = comparison_variants)
    
    # convert into prevalence and frequency estimates
    l[[i]] <- df_match |>
      filter(match_pos) |>
      filter(!ambiguous) |>
      summarise(variant = target_variants[i],
                prev = mean(match),
                freq = mean(prop),
                sample_total = n())
  }
  
  ret <- bind_rows(l)
  return(ret)
}

# ----------------------------------------------------------------------
#' Write a prevalence table to file
#'
#' @description
#' Given a data.frame of prevalence and frequency estimates, writes to a TSV file.
#'
#' @param df_prev a data.frame of prevalence and frequency estimates. Must have
#'   columns variant, prev, freq, sample_total.
#' @param output_path the path to write to.

write_prev <- function(df_prev,
                       output_path) {
  
  # Check input arguments
  stopifnot(is.data.frame(df_prev))
  stopifnot(all(names(df_prev) == c("variant", "prev", "freq", "sample_total")))
  stopifnot(is.character(output_path))
  
  write_tsv(df_prev, file = output_path)
}


# ----------------------------------------------------------------------
# RUN MODULE

# parse arguments
args <- parse_args(OptionParser(option_list = opts))

# read in data and convert to variant string format
variant_strings <- aa_table_to_variant(input_path = args$input_path)

# get the component genotypes that make up these samples
component_variants <- variantstring::get_component_variants(x = variant_strings) |>
  unlist() |>
  na.omit() |>
  as.vector() |>
  unique()

# calculate prevalence
df_prev <- calculate_variant_prevalence(target_variants = component_variants,
                                        comparison_variants = variant_strings)

# write table to file
write_prev(df_prev = df_prev,
           output_path = args$output_path)
