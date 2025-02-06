library(dplyr)
library(optparse)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(validate)

# install a specific version of variantstring package (this is in development so
# may not always be backwards-compatible)
variantstring_version <- "1.8.0"
if (!requireNamespace("variantstring", quietly = TRUE) ||
      packageVersion("variantstring") != variantstring_version) {
   stop(
      "This script requires variantstring version ",  
      variantstring_version, 
      " and it is not installed", 
      call. = FALSE
    )
}
library(variantstring) 

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--aa_table", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing amino acid calls, with the columns: specimen_id, ", 
      "gene, pos, read_count, aa"
    )
  ), 
  make_option(
    "--loci_groups_input", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing loci group definitions, with the ", 
      "columns: group_id, gene_id, aa_position"
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

#' Read amino acid calls table into a tibble
#'
#' @description
#' Takes the path of a TSV file of amino acid calls. Reads these in, 
#' computes the number of amino acids at each position, and returns a 
#' tibble with this information.
#'
#' @param aa_table Path to a TSV file with amino acid calls. It should have
#'   columns for specimen_id, gene, pos, read_count, aa.
#' 
#' @import dplyr
#' 
#' @return Tibble of amino acid calls with specimen_id, gene, pos, 
#'   read_count, aa, and n_aa columns.
create_aa_table_input <- function(aa_table) {
  
  # Check input arguments
  stopifnot(is.character(aa_table))
  
  # read in amino acid calls and validate columns
  df_aa <- read.table(aa_table, header = TRUE)
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
  df_aa <- df_aa |>
    group_by(specimen_id, gene, pos) |>
    mutate(n_aa = n()) |>
    ungroup()

  return(df_aa)
}

#' Read in loci groups table
#'
#' This function takes the path to a table of loci groups for 
#' multilocus allele frequency and prevalence calculations and reads it 
#' into a tibble.
#'
#' @param loci_groups_path Path to loci groups TSV. It should have 
#'   columns for group_id, gene_id, and aa_position.
#'
#' @import dplyr
#'
#' @return Tibble of loci groups, with columns for group_id, gene_id, 
#'   and aa_position.
create_loci_group_input <- function(loci_groups_path) {

  # Check input arguments
  stopifnot(is.character(loci_groups_path))
  
  # Read and validate table
  loci_groups <- read_tsv(
    loci_groups_path, 
    col_types = cols(.default = col_character(), aa_position = col_integer()), 
    progress = FALSE
  )
  rules <- validate::validator(
    is.character(group_id), 
    is.character(gene_id), 
    is.integer(aa_position), 
    ! is.na(group_id), 
    ! is.na(gene_id), 
    ! is.na(aa_position)
  )
  fails <- validate::confront(loci_groups, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  return(loci_groups)

}

#' Convert long form amino acid table to variant strings
#'
#' @description
#' Takes a tibble of amino acid calls. Converts them to 
#' \href{https://github.com/mrc-ide/variantstring}{variant strings},
#' which are then returned.
#'
#' @param aa_table Tibble of amino acid calls. It should have columns 
#'   for specimen_id, gene, pos, read_count, aa, and n_aa.
#' 
#' @import dplyr
#' 
#' @return vector of variant strings corresponding to each specimen.
aa_table_to_variant <- function(aa_table) {
  # get list of long format data.frames
  long_list <- aa_table |>
    mutate(het = (n_aa > 1),
           phased = FALSE) |>
    select(gene, pos, n_aa, het, phased, aa, read_count) |>
    split(f = aa_table$specimen_id)
  names(long_list) <- NULL
  
  # convert to variant strings
  vs <- long_to_variant(long_list)
  
  return(vs)
}


#' Get all component genotypes from a vector of variant strings
#'
#' @description
#' Takes a vector of variant strings and extracts all multi-locus genotypes that
#' can be unambiguously identified within these samples. Return these in variant
#' string format, which will contain no heterozygous loci.
#'
#' @param variant_strings a vector of variant strings.
#' 
#' @return vector of component variant strings.
extract_component_variants <- function(variant_strings) {
  
  # Check input arguments
  check_variant_string(variant_strings)
  
  # get all unique positions
  pos_unique <- position_from_variant_string(variant_strings) |>
    unique()
  
  # get all variant strings at all positions of interest
  vs_allpos <- mapply(function(x) {
    subset_position(position_string = x,
                    variant_strings = variant_strings) |>
      na.omit() |>
      unique()
  }, pos_unique, SIMPLIFY = FALSE) |>
    unlist() |>
    unique()
  
  # get targets as all component genotypes
  targets <- get_component_variants(vs_allpos) |>
    unlist() |>
    na.omit() |>
    unique()
  
  return(targets)
}

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

#' Compute the frequency and prevalence for a group of loci
#'
#' Given one tibble defining the loci in a group of interest and 
#' another containing amino acid calls, this function finds the amino 
#' acid calls corresponding to the group of interest, uses 
#' `aa_table_to_variant()` and `extract_component_variants()` to 
#' transform the data into the formats needed by variantstring, and 
#' then calls `calculate_variant_prevalence()` to compute frequency and 
#' prevalence.
#'
#' @param group_loci A tibble defining a loci group of interest, with 
#'   columns for gene_id and aa_position.
#' @inheritParams aa_table_to_variant
#'
#' @return A tibble with variant, prev, freq, and sample_total columns
compute_prevfreq_for_group <- function(group_loci, aa_table) {
  variant_strings <- aa_table |>
    # Use inner_join() to remove amino acid calls not relevant to loci 
    # in this group
    inner_join(group_loci, by = c("gene" = "gene_id", "pos" = "aa_position")) |>
    aa_table_to_variant()
  component_variants <- extract_component_variants(variant_strings)
  group_prev <- calculate_variant_prevalence(
    component_variants, 
    variant_strings
  )

  return(group_prev)
}

#' Write a prevalence table to file
#'
#' @description
#' Given a data.frame of prevalence and frequency estimates, writes to a TSV file.
#'
#' @param df_prev a data.frame of prevalence and frequency estimates. Must have
#'   columns group_id, variant, prev, freq, sample_total.
#' @param output_path the path to write to.
write_prev <- function(df_prev,
                       output_path) {
  
  # Check input arguments
  stopifnot(is.data.frame(df_prev))
  stopifnot(
    all(
      names(df_prev) == c("group_id", "variant", "prev", "freq", "sample_total")
    )
  )
  stopifnot(is.character(output_path))
  
  write_tsv(df_prev, file = output_path)
}


# RUN MODULE -----------------------------------------------------------

# parse arguments
args <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  args$aa_table <- "../../data/example_amino_acid_calls.tsv"
  args$loci_groups_input <- "../../data/example_loci_groups.tsv"
  args$output_path <- "../../mlafp.tsv"
}


# Read in data
aa_table <- create_aa_table_input(args$aa_table)
loci_groups <- create_loci_group_input(args$loci_groups_input)


# Compute frequency and prevalence for each group
df_prev <- loci_groups |>
  nest(group_loci = c(gene_id, aa_position)) |>
  mutate(group_mlafp = map(group_loci, compute_prevfreq_for_group, aa_table)) |>
  select(-group_loci) |>
  unnest(group_mlafp)

# write table to file
write_prev(df_prev = df_prev, output_path = args$output_path)
