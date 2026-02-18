# Estimate single-locus allele frequency naively with Dcifer

# Load required libraries ----------------------------------------------
library(dcifer)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(dplyr)
library(magrittr)
library(optparse)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing alleles, with the columns: specimen_id, target_id, ", 
      "read_count, and seq. Required."
    )
  ), 
  make_option(
    "--coi_table", 
    help = 
      str_c(
        "TSV containing specimen COIs, with the columns: specimen_id and ", 
        "coi. Optional."
      )
  ), 
  make_option(
    "--slaf_output", 
    help = str_c(
      "Path of TSV file to contain single locus allele frequencies, with the ", 
      "columns: target_id, seq, and freq. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../results/full/mh/mh_pgecore.tsv"
  arg$coi_table <- "../../results/full/coi_table.tsv"
  arg$slaf_output <- "../../results/full/slaf.tsv"
}

#' Read allele table into a tibble
#'
#' Read the allele table TSV into a tibble
#'
#' @param allele_table_path Path to the TSV file containing the allele 
#'   table. There should be character columns for specimen_id, 
#'   target_id, and seq.
#'
#' @return A tibble containing columns for specimen_id, target_id, and 
#'   seq.
create_allele_table_input <- function(allele_table_path) {

  # Read in table
  allele_table <- read_tsv(
    allele_table_path, 
    col_types = cols(.default = col_character()), 
    progress = FALSE
  )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(target_id), 
    is.character(seq), 
    ! is.na(specimen_id), 
    ! is.na(target_id), 
    ! is.na(seq)
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
#'   have a character specimen_id column and an integer coi column.
#' @param allele_list The list format output by `dcifer::readDat` and 
#'   `dcifer::formatDat`.
#'
#' @return Vector of COI values, one for each sample.
create_coi_input <- function(coi_path, allele_list) {

  # Read input table
  coi <- read_tsv(
    coi_path, 
    col_types = cols(.default = col_character(), coi = col_integer()), 
    progress = FALSE
  )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_id), 
    is.integer(coi), 
    ! is.na(specimen_id), 
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
  specimen_inallele_notincoi <- setdiff(names(allele_list), coi$specimen_id)
  if (length(specimen_inallele_notincoi) > 0) {
    stop(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_inallele_notincoi, collapse = " "), 
      call. = FALSE
    )
  }
  specimen_incoi_notinallele <- setdiff(coi$specimen_id, names(allele_list))
  if (length(specimen_incoi_notinallele) > 0) {
    warning(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_incoi_notinallele, collapse = " "), 
      call. = FALSE
    )
  }

  # Join COI to allele list to ensure order matches
  coi <- tibble(specimen_id = names(allele_list)) %>%
    left_join(coi, by = "specimen_id")

  # Create and return named vector of COI
  setNames(coi$coi, coi$specimen_id)

}

#' Convert allele frequency list to tibble and write to TSV
#'
#' This function takes allele frequencies in the list format output by 
#' `dcifer::calcAfreq()`, converts them into a tibble, and writes that 
#' tibble to a TSV.
#'
#' @param allele_freqs_list List of allele frequencies produced by 
#'   `dcifer::calcAfreq()`.
#' @param n_samp_per_target Tibble containing sample sizes for each 
#'   locus, with target_id and sample_total columns.
#' @param output_path Path for TSV of allele frequencies.
write_slaf_output <- function(
                              allele_freqs_list, 
                              n_samp_per_target, 
                              output_path) {
  onetargetaf_list2tib <- function(onetargetaf) {
    tibble(seq = names(onetargetaf), freq = unname(onetargetaf))
  }
  tibble(
      target_id = names(allele_freqs_list), 
      alleles_freqs = unname(allele_freqs_list)
    ) %>%
    mutate(alleles_freqs = map(alleles_freqs, onetargetaf_list2tib)) %>%
    unnest(alleles_freqs) %>%
    left_join(n_samp_per_target, by = "target_id") %>%
    write_tsv(output_path)
}

# Read in allele table -------------------------------------------------
allele_table <- create_allele_table_input(arg$allele_table)
# Convert to Dcifer format and return
dcifer_alleles <- dcifer::formatDat(
    allele_table, 
    svar = "specimen_id", 
    lvar = "target_id", 
    avar = "seq"
  )

# If no COI input was provided, use Dcifer's built-in naive estimation -
if (is.null(arg$coi_table)) {
  coi <- dcifer::getCOI(dcifer_alleles)
} else {
  coi <- create_coi_input(arg$coi_table, dcifer_alleles)
}

# Compute allele frequencies -------------------------------------------
allele_freqs_list <- dcifer::calcAfreq(dcifer_alleles, coi, tol = 1e-5)

# Compute sample sizes for each locus ----------------------------------
n_samp_per_target <- allele_table %>%
  group_by(target_id) %>%
  summarize(sample_total = n_distinct(specimen_id), .groups = "drop")

# Reformat to table and write to disk ----------------------------------
write_slaf_output(allele_freqs_list, n_samp_per_target, arg$slaf_output)
