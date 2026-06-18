#!/usr/bin/env Rscript

# Estimate single-locus allele frequency naively with Dcifer

# Load required libraries ----------------------------------------------
library(dcifer)
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
    "--coi_table", 
    help = 
      str_c(
        "TSV containing specimen COIs, with a coi column and a column with ", 
        "specimen IDs named according to --specimen_name_col. Optional."
      )
  ), 
  make_option(
    "--specimen_name_col", 
    default = "specimen_name", 
    help = "String giving the name of the specimen ID column"
  ), 
  make_option(
    "--target_name_col", 
    default = "target_name", 
    help = 
      "String giving the name of the target name column (e.g., the locus name)"
  ), 
  make_option(
    "--target_value_col", 
    default = "seq", 
    help = 
      str_c(
        "String giving the name of the target value column (e.g., the allele ", 
        "call)"
      )
  ), 
  make_option(
    "--tol", 
    type = "double", 
    default = 1e-04, 
    help = 
      str_c(
        "Double specifying the convergence tolerance. Passed to ", 
        "dcifer::calcAfreq(). Optional."
      )
  ), 
  make_option(
    "--qstart", 
    type = "double", 
    default = 0.5, 
    help = 
      str_c(
        "Double specifying the starting value for frequency estimation. ", 
        "Passed to dcifer::calcAfreq(). Optional."
      )
  ), 
  make_option(
    "--coi_lrank", 
    type = "integer", 
    default = 2, 
    help = 
      str_c(
        "Integer specifying the rank of the locus used to determine ", 
        "sample COI. Passed to dcifer::getCOI(). Should not be combined with ", 
        "--coi_table. Optional."
      )
  ), 
  make_option(
    "--slaf_output", 
    help = str_c(
      "Path of TSV file to contain single locus allele frequencies, with ", 
      "columns with names matching --target_name_col and --target_value_col, ", 
      "and a freq column. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example_allele_table.tsv"
  arg$slaf_output <- "../../slaf.tsv"
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

  # Read in table. Force the specimen name column to character so that
  # all-numeric specimen names are not inferred as numeric (which would
  # drop leading zeros and break downstream string operations).
  allele_table <- read_tsv(
      allele_table_path,
      col_types = do.call(
        cols,
        c(
          setNames(list(col_character()), specimen_name_col),
          list(.default = col_character())
        )
      ),
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

  allele_table

}

#' Read COI table into format needed by Dcifer
#'
#' Read the TSV of specimen COIs into a tibble, join with specimen IDs 
#' from the allele list to ensure order matches, and return as a vector.
#'
#' @param coi_path Path to TSV containing specimen COIs. It should 
#'   have a character specimen_name column and an integer coi column.
#' @param allele_list The list format output by `dcifer::readDat` and 
#'   `dcifer::formatDat`.
#' @inheritParams create_allele_table_input
#'
#' @return Vector of COI values, one for each sample.
create_coi_input <- function(
                             coi_path, 
                             allele_list, 
                             specimen_name_col = "specimen_name") {

  # Read input table. Force the specimen name column to character so that
  # all-numeric specimen names are not inferred as numeric.
  coi <- read_tsv(
      coi_path,
      col_types = do.call(
        cols,
        c(
          setNames(list(col_character()), specimen_name_col),
          list(.default = col_character(), coi = col_integer())
        )
      ),
      progress = FALSE
    ) %>%
    select(all_of(specimen_name_col), coi) %>%
    rename(specimen_name = all_of(specimen_name_col))

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_name), 
    is.integer(coi), 
    ! is.na(specimen_name), 
    ! is.na(coi)
  )
  fails <- validate::confront(coi, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  # Check that all specimen IDs in the allele table are in the COI 
  # input, and vice versa
  specimen_inallele_notincoi <- setdiff(names(allele_list), coi$specimen_name)
  if (length(specimen_inallele_notincoi) > 0) {
    stop(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_inallele_notincoi, collapse = " "), 
      call. = FALSE
    )
  }
  specimen_incoi_notinallele <- setdiff(coi$specimen_name, names(allele_list))
  if (length(specimen_incoi_notinallele) > 0) {
    warning(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_incoi_notinallele, collapse = " "), 
      call. = FALSE
    )
  }

  # Join COI to allele list to ensure order matches
  coi <- tibble(specimen_name = names(allele_list)) %>%
    left_join(coi, by = "specimen_name")

  # Create and return named vector of COI
  setNames(coi$coi, coi$specimen_name)

}

#' Convert allele frequency list to tibble
#'
#' This function takes allele frequencies in the list format output by 
#' `dcifer::calcAfreq()` and converts them into a tibble in preparation 
#' for writing to disk.
#'
#' @param allele_freqs_list List of allele frequencies produced by 
#'   `dcifer::calcAfreq()`.
#' @param n_samp_per_target Tibble containing sample sizes for each 
#'   locus, with target_name and sample_total columns.
#' @inheritParams create_allele_table_input
prepare_slaf_output <- function(
                                allele_freqs_list, 
                                n_samp_per_target, 
                                target_name_col = "target_name", 
                                target_value_col = "seq") {
  onetargetaf_list2tib <- function(onetargetaf) {
    tibble(target_value = names(onetargetaf), freq = unname(onetargetaf))
  }
  tibble(
      target_name = names(allele_freqs_list), 
      alleles_freqs = unname(allele_freqs_list)
    ) %>%
    mutate(alleles_freqs = map(alleles_freqs, onetargetaf_list2tib)) %>%
    unnest(alleles_freqs) %>%
    left_join(n_samp_per_target, by = "target_name") %>%
    # Revert column names back to user-specified names
    rename(
      !!target_name_col := target_name, 
      !!target_value_col := target_value
    )
}

# Read in allele table -------------------------------------------------
allele_table <- create_allele_table_input(
  arg$allele_table, 
  specimen_name_col = arg$specimen_name_col, 
  target_name_col = arg$target_name_col, 
  target_value_col = arg$target_value_col
)
# Convert to Dcifer format and return
dcifer_alleles <- dcifer::formatDat(
    allele_table, 
    svar = "specimen_name", 
    lvar = "target_name", 
    avar = "target_value"
  )

# If no COI input was provided, use Dcifer's built-in naive estimation -
if (is.null(arg$coi_table)) {
  coi <- dcifer::getCOI(dcifer_alleles, lrank = arg$coi_lrank)
} else {
  if (arg$coi_lrank != 2) {
    warning(
      "The --coi_lrank argument has been changed from the default, but this ", 
      "will have no effect as --coi_table has also been specified", 
      call. = FALSE
    )
  }
  coi <- create_coi_input(
    arg$coi_table, 
    dcifer_alleles, 
    specimen_name_col = arg$specimen_name_col
  )
}

# Compute allele frequencies -------------------------------------------
allele_freqs_list <- dcifer::calcAfreq(
    dcifer_alleles, 
    coi, 
    tol = arg$tol, 
    qstart = arg$qstart
  )

# Compute sample sizes for each locus ----------------------------------
n_samp_per_target <- allele_table %>%
  group_by(target_name) %>%
  summarize(sample_total = n_distinct(specimen_name), .groups = "drop")

# Reformat to table and write to disk ----------------------------------
prepare_slaf_output(
    allele_freqs_list, 
    n_samp_per_target, 
    target_name_col = arg$target_name_col, 
    target_value_col = arg$target_value_col
 ) %>%
  write_tsv(arg$slaf_output)
