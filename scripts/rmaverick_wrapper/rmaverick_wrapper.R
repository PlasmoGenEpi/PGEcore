#!/usr/bin/env Rscript

# Estimate genetic clusters with rmaverick

# Load required libraries ----------------------------------------------
library(rmaverick)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(optparse)
library(parallel)
library(readr)
library(stringr)
library(tidyr, warn.conflicts = FALSE)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing alleles, with columns identifying specimens, ", 
      "target names, and target values. The names of these columns are given ", 
      "by the --specimen_name_col, --target_name_col, and --target_value_col ", 
      "arguments, respectively. Required."
    )
  ), 
  make_option(
    "--specimen_name_col", 
    default = "specimen_name", 
    help = "String giving the name of the specimen ID column. Optional."
  ), 
  make_option(
    "--target_name_col", 
    default = "target_name", 
    help = 
      str_c(
        "String giving the name of the target name column (e.g., the locus ", 
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
    "--Kmax", 
    default = 5, 
    type = "integer", 
    help = "Largest K value to evaluate. Optional."
  ), 
  make_option(
    "--burnin", 
    default = 1000, 
    type = "integer", 
    help = "The number of burn-in iterations. Optional."
  ), 
  make_option(
    "--samples", 
    default = 1000, 
    type = "integer", 
    help = "The number of sampling iterations. Optional."
  ), 
  make_option(
    "--rungs", 
    default = 10, 
    type = "integer", 
    help = "The number of temperature rungs. Optional."
  ), 
  make_option(
    "--GTI_pow", 
    default = 3, 
    type = "double", 
    help = str_c(
      "The power used in the generalised thermodynamic integration method. ", 
      "Must be greater than 1.1. Optional."
    )
  ), 
  make_option(
    "--coupling_on", 
    action = "store_true", 
    default = FALSE, 
    help = str_c(
      "Whether to implement Metropolis-coupling over temperature rungs. ",
      "Optional."
    )
  ), 
  make_option(
    "--threads", 
    default = 1, 
    type = "integer", 
    help = "Number of threads to use. Optional."
  ), 
  make_option(
    "--seed", 
    default = 1, 
    type = "integer", 
    help = "Seed for random number generation. Optional."
  ), 
  make_option(
    "--model_results_output", 
    help = str_c(
      "Path of RDS file to contain the rmaverick results object. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table_monos_only.tsv"
  arg$model_results_output <- "../../rmaverick_res.rds"
  arg$threads <- 5
  arg$Kmax <- 5
}

#' Read allele table into a tibble
#'
#' Read the allele table TSV into a tibble
#'
#' @param allele_table_path Path to the TSV file containing the allele 
#'   table.
#' @param specimen_name_col String giving the name of the specimen ID 
#'   column.
#' @param target_name_col String giving the name of the target name column.
#' @param target_value_col String giving the name of the column 
#'   containing target values (i.e., the genotypes).
#'
#' @return A tibble containing columns for specimen_name, target_name, and 
#'   target_value.
create_allele_table_input <- function(
                                      allele_table_path, 
                                      specimen_name_col = "specimen_name", 
                                      target_name_col = "target_name", 
                                      target_value_col = "seq") {

  # Read in table
  allele_table <- read_tsv(
      allele_table_path, 
      col_types = cols(.default = col_character()), 
      progress = FALSE
    ) %>%
    select(all_of(c(specimen_name_col, target_name_col, target_value_col))) %>%
    # Standardize names
    rename(
      specimen_name = all_of(specimen_name_col), 
      target_name = all_of(target_name_col), 
      target_value = all_of(target_value_col)
    )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_name), 
    is.character(target_name), 
    is.character(target_value), 
    ! is.na(specimen_name), 
    ! is.na(target_name), 
    ! is.na(target_value)
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

  # Filter out invariant loci
  allele_table <- allele_table %>%
    group_by(target_name) %>%
    filter(n_distinct(target_value) > 1)

  # Recode microhaplotypes as arbitrary integers
  allele_table <- allele_table %>%
    group_by(target_name) %>%
    mutate(haplotype = as.integer(factor(target_value))) %>%
    ungroup()
  # Store map for recovery later
  target_value_map <- allele_table %>%
    distinct(target_name, target_value, haplotype)
  allele_table <- allele_table %>%
    select(-target_value)

  # Make implicit missing data explicit and code as -9
  allele_table <- allele_table %>%
    complete(specimen_name, target_name, fill = list(haplotype = -9))

  # Pivot to wide format
  allele_table <- allele_table %>%
    pivot_wider(
      names_from = target_name, 
      values_from = haplotype
    ) %>%
    # Add arbitary population and ploidy columns
    mutate(population = 1, ploidy = 1) %>%
    relocate(population, ploidy, .after = 1) %>%
    # Necessary to avoid cryptic bugs from rmaverick
    as.data.frame()

  return(
    list(
      allele_table = allele_table, 
      target_value_map = target_value_map
    )
  )

}

#' Create rmaverick project and run MCMC
#'
#' This function takes an allele table of the format expected by 
#' `rmaverick::bind_data()`, creates an rmaverick project with 
#' this data, runs the MCMC to fit the model, and returns the rmaverick 
#' project object.
#'
#' @param allele_data An allele table of the format expected by 
#'   `rmaverick::bind_data()`.
#' @param Kmax Numeric specifying the largest K to evaluate.
#' @param threads Number of threads to use when fitting models for 
#'   multiple K values.
#' @param ... Arguments passed on to `rmaverick::run_mcmc()`.
#'
#' @return An `rmaverick::mavproject()` containing the results.
run_rmaverick <- function(
                          allele_data, 
                          Kmax = 5, 
                          threads = Kmax, 
                          ...) {

  # Set up project
  n_samp <- nrow(allele_data)
  mavproj <- rmaverick::mavproject() %>%
    rmaverick::bind_data(
      df = allele_data, 
      ID_col = 1, 
      pop_col = 2, 
      ploidy_col = 3
    ) %>%
    rmaverick::new_set(
      name = "rmaverick results", 
      admix_on = TRUE
    )

  # Run MCMC, with one K value per thread
  if (threads > 1) {
    cl <- makeCluster(threads)
  } else {
    cl <- NULL
  }
  mavproj <- mavproj %>%
    rmaverick::run_mcmc(
      K = 1:Kmax, 
      cluster = cl, 
      ...
    )
  if (threads > 1) {
    stopCluster(cl)
  }

  mavproj
}

set.seed(arg$seed)

# Read in allele table -------------------------------------------------
alleles_and_maps <- create_allele_table_input(
  arg$allele_table, 
  specimen_name_col = arg$specimen_name_col, 
  target_name_col = arg$target_name_col, 
  target_value_col = arg$target_value_col
)

# Run rmaverick --------------------------------------------------------
rmaverick_res <- run_rmaverick(
  alleles_and_maps$allele_table, 
  Kmax = arg$Kmax, 
  threads = arg$threads, 
  burnin = arg$burnin, 
  samples = arg$samples, 
  rungs = arg$rungs, 
  GTI_pow = arg$GTI_pow, 
  coupling_on = arg$coupling_on
)

# Save results ---------------------------------------------------------
rmaverick_res %>%
  write_rds(arg$model_results_output)
