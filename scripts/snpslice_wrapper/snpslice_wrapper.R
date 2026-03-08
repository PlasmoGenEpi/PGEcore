#!/usr/bin/env Rscript

# Estimate multi-locus allele frequency and COI with SNP-Slice

# Load required libraries ----------------------------------------------
library(snp.slicer)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr, warn.conflicts = FALSE)
library(magrittr)
library(optparse)
library(purrr, warn.conflicts = FALSE)
library(readr)
library(stringr)
library(tibble)
library(tidyr, warn.conflicts = FALSE)

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
    default = "aa_locus", 
    help = 
      str_c(
        "String giving the name of the target ID column (e.g., the locus ", 
        "name). Optional."
      )
  ), 
  make_option(
    "--target_value_col", 
    default = "aa", 
    help = 
      str_c(
        "String giving the name of the target value column (e.g., the allele ", 
        "call). Optional."
      )
  ), 
  make_option(
    "--target_count_col", 
    default = "read_count", 
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
    "--loci_limit", 
    type = "integer", 
    help = str_c(
      "To prevent SNP-Slice from running for a long time on large datasets, ", 
      "the input loci can be subsetted to contain the loci of interest in ", 
      "loci_groups_input and additional random loci, up to the limit ", 
      "specified here. Using more loci will improve COI estimation, but ", 
      "increase computation time. If not provided, all loci will be used."
    )
  ), 
  make_option(
    "--model", 
    default = "negative_binomial", 
    help = str_c(
      "Observation model to use. Options: 'poisson', ", 
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
    "--use_mcmc_for_af_and_coi", 
    action = "store_true", 
    default = FALSE, 
    help = str_c(
      "For calculating COI and allele frequencies, whether to sample from ", 
      "MCMC (TRUE) or use MAP estimates (FALSE; default)."
    )
  ), 
  make_option(
    c("-v", "--verbose"), 
    action = "store_true", 
    default = FALSE, 
    help = "Verbose output."
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
parser <- OptionParser(option_list = opts)
arg <- parse_args(parser)
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_amino_acid_calls.tsv"
  arg$loci_groups_input <- "../../data/example_loci_groups.tsv"
  arg$target_id_col <- "aa_locus"
  arg$target_value_col <- "aa"
  arg$target_count_col <- "read_count"
  arg$loci_limit <- 20
  arg$n_mcmc <- 100
  arg$verbose <- TRUE
  arg$use_mcmc_for_af_and_coi <- TRUE
  arg$mlaf_output <- "../../mlaf.tsv"
  arg$coi_output <- "../../coi.tsv"
}

set.seed(arg$seed)

#' Check for required arguments, and report which are missing 
#'
#' @param parser the parser created from optparse
#' @param arg the parsed arguments from optparse
#' @param required_args the required arguments (without the --)
#'
#' @return returns void if all required arguments
checkOptparseRequiredArgsThrow <- function(parser, arg, required_args){
  missing <- setdiff(required_args, names(arg))
  if(length(missing) > 0){
    missing = paste0("--", missing)
    print_help(parser)
    stop(paste0("missing the following arguments: ", paste0(missing, collapse = ", ")))
  }
}

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
#' @param allele_table A tibble produced by 
#'   `create_allele_table_input()`.
#' @inheritParams create_allele_table_input
#'
#' @import dplyr
#'
#' @return Tibble of loci groups, with columns for group_id and 
#'   target_id.
create_loci_group_input <- function(
                                    loci_groups_path, 
                                    allele_table, 
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

  # Convert to list format
  loci_groups <- split(loci_groups$target_id, loci_groups$group_id)

  # Check for missing and non-biallelic loci
  for (lg in names(loci_groups)) {
    missing_trgs <- setdiff(loci_groups[[lg]], allele_table$target_id)
    if (length(missing_trgs) > 0) {
      warning(
        "The target(s) ", 
        str_c(missing_trgs, collapse = ", "), 
        " in the group ", 
        lg, 
        " are missing and will be excluded.", 
        call. = FALSE
      )
    }
    non_biallelic_trgs <- allele_table %>%
      filter(target_id %in% loci_groups[[lg]]) %>%
      group_by(target_id) %>%
      filter(n_distinct(target_value) > 2) %$%
      unique(target_id)
    if (length(non_biallelic_trgs > 0)) {
      warning(
        "The target(s) ", 
        str_c(non_biallelic_trgs, collapse = ", "), 
        " in the group ", 
        lg, 
        " have more than two alleles and will be excluded.", 
        call. = FALSE
      )
    }
    loci_groups[[lg]] <- setdiff(
      loci_groups[[lg]], 
      c(missing_trgs, non_biallelic_trgs)
    )
  }
  if (length(unlist(loci_groups)) == 0) {
    stop("The data is missing all loci in loci groups", call. = FALSE)
  }

  return(loci_groups)

}

#' Calculate and format allele frequencies from SNP-Slice results
#'
#' This function takes a snp.slicer results object and a list of loci 
#' groups, calculates the allele frequencies for each group, and 
#' formats the output into a tibble suitable for writing to disk.
#'
#' @param snp_slice_res A snp.slicer results object produced by 
#'   `snp.slicer::snp_slice()`.
#' @param loci_groups A list containing named character vectors 
#'   defining loci groups.
#' @param use_mcmc Logical indicating whether to use the MCMC results 
#'   for calculating allele frequencies.
#'
#' @return A tibble with group_id, variant, freq, allele_total, and 
#'   sample_total columns.
prepare_af_output <- function(
                              snp_slice_res, 
                              loci_groups, 
                              use_mcmc) {

  # Reformat a tibble of allele frequencies to have variant string names
  format_af_table_w_variantstring <- function(af_table, group_id, loci_groups) {
    # Format a multi-locus genotype into long form for variantstring
    prep_variantstring_input <- function(allele, loci_names) {
      aa <- str_split_1(allele, "\\|")
      tibble(gene_pos = loci_names, aa = aa) %>%
        # SNP-Slice outputs frequencies for some haplotypes that 
        # include ?s to indicate unknown alleles at that position. 
        # There is no clear way to format these as a variant string, so 
        # instead these positions are omitted and the smaller, 
        # unambiguous haplotype is reported.
        filter(aa != "?") %>%
        separate_wider_delim(gene_pos, ":", names = c("gene", "pos")) %>%
        mutate(pos = as.integer(pos)) %>%
        mutate(n_aa = 1, het = FALSE, phased = TRUE, read_count = NA) %>%
        relocate(aa, .before = read_count)
    }

    # Convert allele names to variant string format
    af_table %>%
      as_tibble() %>%
      # All possible genotypes will be included, but many will have a 
      # freq of 0 if use_mcmc = FALSE or some samples are missing loci
      filter(frequency > 0) %>%
      mutate(
        allele = map(
          allele, 
          prep_variantstring_input, 
          loci_groups[[group_id]]
        )
      ) %>%
      # After the filtering in prep_variantstring_input(), it is 
      # possible allele will be an empty string if all of the positions 
      # were unknown (i.e., a "?")
      filter(map_lgl(allele, ~nrow(.x) > 0)) %>%
      mutate(allele = variantstring::long_to_variant(allele))
  }

  # Calculate allele frequencies
  snp_slicer_af_by_group <- snp.slicer::calculate_allele_frequencies_by_sets(
      snpslice_res, 
      loci_groups, 
      use_map = ! use_mcmc
    )
  # Reformat
  af_tib <- tibble(
      group_id = names(snp_slicer_af_by_group), 
      af_tib = snp_slicer_af_by_group
    ) %>%
    mutate(
      af_tib = map2(
        af_tib, 
        group_id, 
        format_af_table_w_variantstring, 
        loci_groups
      )
    ) %>%
    unnest(af_tib) %>%
    rename(
      variant = allele, 
      freq = frequency
    )
  if (use_mcmc) {
    af_tib <- af_tib %>%
      select(-mean_count, -n_samples)
  } else {
    af_tib <- af_tib %>%
      rename(allele_total = total_parasites, allele_count = count)
  }

  return(af_tib)

}

#' Calculate and format COI from SNP-Slice results
#'
#' This function takes a snp.slicer results object, calculates COI from 
#' this object, and formats the output into a tibble suitable for 
#' writing to disk.
#'
#' @param snp_slice_res A snp.slicer results object produced by 
#'   `snp.slicer::snp_slice()`.
#' @inheritParams create_allele_table_input
#' @inheritParams prepare_af_output
#'
#' @return A tibble with a coi column and a specimen ID column with name 
#'   matching specimen_id_col.
prepare_coi_output <- function(
                              snp_slice_res, 
                              specimen_id_col, 
                              use_mcmc) {
  coi_tib <- snp_slice_res %>%
    snp.slicer::calculate_individual_coi(
      use_map = ! arg$use_mcmc_for_af_and_coi
    ) %>%
    select(-host_index) %>%
    rename(!! arg$specimen_id_col := host_id, coi = coi_estimate)
  if (! arg$use_mcmc_for_af_and_coi) {
    coi_tib <- coi_tib %>%
      select(-coi_sd, -coi_lower, -coi_upper)
  }

  return(coi_tib)
}

# Check for required arguments -----------------------------------------
required_arguments = c(
  "allele_table", 
  "loci_groups_input", 
  "mlaf_output", 
  "coi_output"
)
checkOptparseRequiredArgsThrow(parser, arg, required_arguments)

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
  allele_table, 
  target_id_col = arg$target_id_col
)

# Subset loci if loci_limit provided -----------------------------------
if (! is.null(arg$loci_limit)) {
  if (n_distinct(allele_table$target_id) > arg$loci_limit) {
    # Loci in loci groups that must be included
    loci_oi <- unique(unlist(loci_groups))
    # Randomly sample additional loci to reach limit
    n_loci_select <- arg$loci_limit - length(loci_oi)
    loci_random <- sample(
      setdiff(allele_table$target_id, loci_oi), 
      n_loci_select
    )
    loci_selected <- c(loci_oi, loci_random)
    allele_table <- allele_table %>%
      filter(target_id %in% loci_selected)
  }
}

# Run SNP-Slice --------------------------------------------------------
snpslice_res <- snp.slicer::snp_slice(
  allele_table, 
  model = arg$model, 
  n_mcmc = arg$n_mcmc, 
  burnin = arg$burnin, 
  alpha = arg$alpha, 
  rho = arg$rho, 
  threshold = arg$threshold, 
  gap = arg$gap, 
  store_mcmc = TRUE, 
  verbose = arg$v, 
  specimen_id_col = "specimen_id", 
  target_id_col = "target_id", 
  target_value_col = "target_value", 
  target_count_col = "target_count"
)

# Calculate and write allele frequencies -------------------------------
snpslice_res %>%
  prepare_af_output(loci_groups, arg$use_mcmc_for_af_and_coi) %>%
  write_tsv(arg$mlaf_output)

# Calculate and write COI ----------------------------------------------
snpslice_res %>%
  prepare_coi_output(arg$specimen_id_col, arg$use_mcmc_for_af_and_coi) %>%
  write_tsv(arg$coi_output)
