#!/usr/bin/env Rscript

# Estimate COI, allele frequencies, and genetic clustering with MALECOT

# Load required libraries ----------------------------------------------
library(MALECOT)
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
    "--COI_model", 
    default = "nb", 
    help = str_c(
      'the type of prior on COI. Must be one of "uniform", "poisson", or ', 
      '"nb" (negative binomial). Optional.'
    )
  ), 
  make_option(
    "--COI_max", 
    default = 20, 
    type = "integer", 
    help = "The maximum COI allowed for any given sample. Optional."
  ), 
  make_option(
    "--use_provided_mean_COI", 
    action = "store_true", 
    default = FALSE, 
    help = str_c(
      "Whether to use mean COI provided or have MALECOT estimate it. Optional."
    )
  ), 
  make_option(
    "--COI_mean", 
    default = 3, 
    type = "double", 
    help = str_c(
      "Single scalar value specifying the mean COI for all subpopulations. ", 
      "Optional."
    )
  ), 
  make_option(
    "--COI_dispersion", 
    default = 2, 
    type = "double", 
    help = str_c(
      "The ratio of the variance to the mean of the prior on COI. Only ", 
      "applies under the negative binomial model. Must be >1, as a ratio of 1 ", 
      "can be achieved by using the Poisson distribution. Optional."
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
      "Path of RDS file to contain the MALECOT results object. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table.tsv"
  arg$model_results_output <- "../../MALECOT_res.rds"
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
      sample_ID = all_of(specimen_name_col), 
      target_name = all_of(target_name_col), 
      target_value = all_of(target_value_col)
    )

  # Validate fields
  rules <- validate::validator(
    is.character(sample_ID), 
    is.character(target_name), 
    is.character(target_value), 
    ! is.na(sample_ID), 
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

  # Recode locus and haplotype as arbitrary integers
  allele_table <- allele_table %>%
    mutate(locus = as.integer(factor(target_name))) %>%
    group_by(locus) %>%
    mutate(haplotype = as.integer(factor(target_value))) %>%
    ungroup()
  # Store maps for recovery later
  target_name_map <- allele_table %>%
    distinct(target_name, locus)
  target_value_map <- allele_table %>%
    distinct(target_name, target_value, haplotype)
  allele_table <- allele_table %>%
    select(-target_name, -target_value)

  # Make implicit missing data explicit and code as -9
  allele_table <- allele_table %>%
    complete(sample_ID, locus, fill = list(haplotype = -9))

  return(
    list(
      allele_table = allele_table, 
      target_name_map = target_name_map, 
      target_value_map = target_value_map
    )
  )

}

#' Create MALECOT project and run MCMC
#'
#' This function takes an allele table of the format expected by 
#' `MALECOT::bind_data_multiallelic()`, creates a MALECOT project with 
#' this data, runs the MCMC to fit the model, and returns the MALECOT 
#' project object.
#'
#' @param allele_data An allele table of the format expected by 
#'   `MALECOT::bind_data_multiallelic()`.
#' @inheritParams MALECOT::new_set
#' @param Kmax Numeric specifying the largest K to evaluate.
#' @param threads Number of threads to use when fitting models for 
#'   multiple K values.
#' @param ... Arguments passed on to `MALECOT::run_mcmc()`.
#'
#' @return A `MALECOT::malecot_project()` containing the results.
run_malecot <- function(
                        allele_data, 
                        COI_model = "nb", 
                        COI_max = 20, 
                        estimate_COI_mean = TRUE, 
                        COI_mean = 3, 
                        COI_dispersion = 2, 
                        Kmax = 5, 
                        threads = Kmax, 
                        ...) {

  # Set up project
  n_samps <- n_distinct(allele_data$sample_ID)
  malproj <- MALECOT::malecot_project() %>%
    MALECOT::bind_data_multiallelic(df = allele_data) %>%
    MALECOT::new_set(
      name = "MALECOT results", 
      COI_manual = rep(1, n_samps), 
      COI_model = COI_model, 
      COI_max = COI_max, 
      estimate_COI_mean = estimate_COI_mean, 
      COI_mean = COI_mean, 
      COI_dispersion = COI_dispersion, 
      estimate_error = TRUE
    )

  # Run MCMC, with one K value per thread
  if (threads > 1) {
    cl <- makeCluster(threads)
  } else {
    cl <- NULL
  }
  malproj <- malproj %>%
    MALECOT::run_mcmc(
      K = 1:Kmax, 
      cluster = cl, 
      ...
    )
  if (threads > 1) {
    stopCluster(cl)
  }

  malproj
}

set.seed(arg$seed)

# Read in allele table -------------------------------------------------
alleles_and_maps <- create_allele_table_input(
  arg$allele_table, 
  specimen_name_col = arg$specimen_name_col, 
  target_name_col = arg$target_name_col, 
  target_value_col = arg$target_value_col
)

# Run MALECOT ----------------------------------------------------------
malecot_res <- run_malecot(
  alleles_and_maps$allele_table, 
  COI_model = arg$COI_model, 
  COI_max = arg$COI_max, 
  estimate_COI_mean = ! arg$use_provided_mean_COI, 
  COI_mean = arg$COI_mean, 
  COI_dispersion = arg$COI_dispersion, 
  Kmax = arg$Kmax, 
  threads = arg$threads, 
  burnin = arg$burnin, 
  samples = arg$samples, 
  rungs = arg$rungs, 
  GTI_pow = arg$GTI_pow, 
  coupling_on = arg$coupling_on
)

# Save results ---------------------------------------------------------
malecot_res %>%
  write_rds(arg$model_results_output)
