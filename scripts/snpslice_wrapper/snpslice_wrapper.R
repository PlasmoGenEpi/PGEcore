#!/usr/bin/env Rscript

# Estimate multi-locus allele frequency and COI with SNP-Slice

# Load required libraries ----------------------------------------------
library(snp.slicer)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(optparse)
library(readr)
library(stringr)
library(tibble)
library(tidyr, warn.conflicts = FALSE)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing alleles, with columns identifying specimens, ", 
      "target names, target values, and target counts. The names of these ", 
      "columns are given by the --specimen_id_col, --target_id_col, ", 
      "--target_value_col, and --target_count_col arguments, respectively. ", 
      "Required."
    )
  ), 
  make_option(
    "--specimen_id_col", 
    default = "specimen_id", 
    help = "String giving the name of the specimen ID column. Optional."
  ), 
  make_option(
    "--target_id_col", 
    default = "target_id", 
    help = 
      str_c(
        "String giving the name of the target ID column (e.g., the locus ", 
        "name). Optional."
      )
  ), 
  make_option(
    "--target_value_col", 
    default = "seq", 
    help = 
      str_c(
        "String giving the name of the target value column (e.g., the allele ", 
        "call). Optional."
      )
  ), 
  make_option(
    "--target_count_col", 
    default = "seq", 
    help = 
      str_c(
        "String giving the name of the target count column (e.g., the read ", 
        "counts). Optional."
      )
  ), 
  make_option(
    "--loci_groups_input", 
    type = "character",
    help = str_c(
      "Path to a TSV file containing loci group definitions, with a group_id ", 
      "column and a column with name matching --target_id_col. Required."
    )
  ), 
  make_option(
    "--model", 
    default = "negative_binomial", 
    help = str_c(
      "Observation model to use. Options: 'categorical', 'poisson', ", 
      "'binomial', 'negative_binomial' (default). Optional."
    )
  ), 
  make_option(
    "--n_mcmc", 
    type = "integer", 
    default = 10000, 
    help = "Number of MCMC iterations. Optional."
  ), 
  make_option(
    "--burnin", 
    type = "double", 
    default = NULL, 
    help = "Burn-in period. If NULL, defaults to n_mcmc/2. Optional."
  ), 
  make_option(
    "--alpha", 
    type = "double", 
    default = 2.6, 
    help = "IBP concentration parameter. Optional."
  ), 
  make_option(
    "--rho", 
    type = "double", 
    default = 0.5, 
    help = "Dictionary sparsity parameter. Optional."
  ), 
  make_option(
    "--threshold", 
    type = "double", 
    default = 0.001, 
    help = "Threshold for identifying single infections. Optional."
  ), 
  make_option(
    "--gap", 
    type = "integer", 
    default = NULL, 
    help = str_c(
      "Early stopping threshold. If NULL, runs for full n_mcmc iterations. ", 
      "Optional."
    )
  ), 
  make_option(
    "--store_mcmc", 
    action = "store_true", 
    default = FALSE, 
    help = "Whether to store full MCMC samples (default: FALSE)"
  ), 
  make_option(
    "--seed", 
    type = "integer", 
    default = 1, 
    help = "Random number seed. Optional."
  ), 
  make_option(
    "--mlaf_output", 
    help = str_c(
      "Path of TSV file to contain multilocus allele frequencies, with a ", 
      "group_id column, a variant column using the variantstring format, ", 
      "and a freq column. Required."
    )
  ), 
  make_option(
    "--coi_output", 
    help = str_c(
      "Path of TSV file to contain COI estimates, with a specimen_id column ", 
      "and a coi column. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example_amino_acid_calls.tsv"
  arg$loci_groups_input <- "../../data/example_loci_groups.tsv"
  arg$target_id_col <- "aa_locus"
  arg$target_value_col <- "aa"
  arg$target_count_col <- "read_count"
  arg$n_mcmc <- 100
}

#' Read allele table into a tibble
#'
#' Read the allele table TSV into a tibble
#'
#' @param allele_table_path Path to the TSV file containing the allele 
#'   table.
#' @param specimen_id_col String giving the name of the specimen ID 
#'   column.
#' @param target_id_col String giving the name of the target ID column.
#' @param target_value_col String giving the name of the column 
#'   containing target values (i.e., the genotypes).
#' @param target_count_col String giving the name of the column 
#'   containing target counts.
#'
#' @return A tibble containing columns for specimen_id, target_id, 
#'   target_value, and target_count.
create_allele_table_input <- function(
                                      allele_table_path, 
                                      specimen_id_col = "specimen_id", 
                                      target_id_col = "aa_locus", 
                                      target_value_col = "aa", 
                                      target_count_col = "read_count") {

  # Read in table
  allele_table <- read_tsv(
      allele_table_path, 
      col_types = cols(
        .default = col_character(), 
        !!target_count_col := col_double()
      ), 
      progress = FALSE
    ) %>%
    select(
      all_of(
        c(specimen_id_col, target_id_col, target_value_col, target_count_col)
      )
    ) %>%
    # Standardize names
    rename(
      specimen_id = all_of(specimen_id_col), 
      target_id = all_of(target_id_col), 
      target_value = all_of(target_value_col), 
      target_count = all_of(target_count_col)
    )

    print(allele_table)
  # Validate fields
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(target_id), 
    is.character(target_value), 
    is.double(target_count), 
    ! is.na(specimen_id), 
    ! is.na(target_id), 
    ! is.na(target_value), 
    ! is.na(target_count)
  )
  fails <- validate::confront(allele_table, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  allele_table

}

#' Read in loci groups table
#'
#' This function takes the path to a table of loci groups for 
#' multilocus allele frequency and prevalence calculations and reads it 
#' into a tibble.
#'
#' @param loci_groups_path Path to loci groups TSV. It should have 
#'   a column for group_id and a column matching target_id_col.
#' @inheritParams create_allele_table_input
#'
#' @import dplyr
#'
#' @return Tibble of loci groups, with columns for group_id and 
#'   target_id.
create_loci_group_input <- function(
                                    loci_groups_path, 
                                    target_id_col = "target_id") {

  # Check input arguments
  stopifnot(is.character(loci_groups_path))
  
  # Read and validate table
  loci_groups <- read_tsv(
      loci_groups_path, 
      col_types = cols(.default = col_character(), aa_position = col_integer()), 
      progress = FALSE
    ) %>%
    select(group_id, all_of(target_id_col)) %>%
    # Standardize names
    rename(target_id = all_of(target_id_col))
  rules <- validate::validator(
    is.character(group_id), 
    is.character(target_id), 
    ! is.na(group_id), 
    ! is.na(target_id)
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

options(error = stop)
# Read inputs ----------------------------------------------------------
allele_table <- create_allele_table_input(
  arg$allele_table, 
  specimen_id_col = arg$specimen_id_col, 
  target_id_col = arg$target_id_col, 
  target_value_col = arg$target_value_col, 
  target_count_col = arg$target_count_col
)
loci_groups <- create_loci_group_input(
  arg$loci_groups_input, 
  target_id_col = arg$target_id_col
)

snpslice_res <- snp.slicer::snp_slice(
  allele_table, 
  model = arg$model, 
  n_mcmc = arg$n_mcmc, 
  burnin = arg$burnin, 
  alpha = arg$alpha, 
  rho = arg$rho, 
  threshold = arg$threshold, 
  gap = arg$gap, 
  store_mcmc = arg$store_mcmc, 
  specimen_id_col = "specimen_id", 
  target_id_col = "target_id", 
  target_value_col = "target_value", 
  target_count_col = "target_count"
)
