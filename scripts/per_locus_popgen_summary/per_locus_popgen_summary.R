#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(optparse)
library(stringr)

library(ape)
library(msa)
library(pegas)

# Parse arguments ------------------------------------------------------
opts = list(
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
    "--out",
    default = "per_locus_popgen_summary.tsv",
    help = str_c(
      "the output path of the results"
    )
  ), 
  make_option(
    "--msa_method",
    default = "Muscle",
    help = str_c(
      "the default is %default, options are 'ClustalW', 'ClustalOmega', 'Muscle'"
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table.tsv"
  arg$out <- "../../popgen_summary.tsv"
}

# locus counting functions -----------------------------------------------------

# create locus data ------------------------------------------------------------
#' Create Locus Data
#'
#' Reads an input file, validates its format, and creates a locus data frame
#' containing sample IDs, target IDs, and allele sequences.
#'
#' @param input_path A string specifying the path to the input file.
#' @param specimen_name_col String giving the name of the specimen ID 
#'   column.
#' @param target_name_col String giving the name of the target name column.
#' @param target_value_col String giving the name of the column 
#'   containing target values (i.e., the genotypes).
#'
#' @return A tibble containing columns for specimen_name, target_name, and 
#'   target_value.
create_locus_data <- function(
                              input_path, 
                              specimen_name_col = "specimen_name", 
                              target_name_col = "target_name", 
                              target_value_col = "seq") {

  input_data <- read_tsv(
      input_path, 
      col_types = cols(
        .default = col_character(), 
        reads = col_double()
      ), 
      progress = FALSE
    )
  locus_data <- input_data %>%
    select(all_of(c(specimen_name_col, target_name_col, target_value_col))) %>%
    # Standardize names
    rename_with(
      ~ c("specimen_name", "target_name", "target_value"),
      .cols = all_of(c(specimen_name_col, target_name_col, target_value_col))
    )
  
  rules <- validate::validator(
    # Data columns
    is.character(specimen_name),
    is.character(target_name),
    is.character(target_value),
  
    # Non-missing values
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(target_value)
  )
  
  # Confront the analysis_object with validation rules
  fails <- validate::confront(locus_data, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  
  # Raise an error if any validations fail
  if (nrow(fails) > 0) {
    stop(
      "Analysis object failed one or more validation checks: ",
      str_c(fails$expression, collapse = "\n"),
      call. = FALSE
    )
  }
  
  return(locus_data)
}

# calculate popgen stats -------------------------------------------------------
# Calculate popgen stats 
#' Calculate Population Genetics Statistics
#'
#' Computes population genetic statistics, including nucleotide diversity,
#' the number of segregating sites, and Tajima's D, for a set of allele sequences.
#'
#' @param allele_data A character vector containing allele sequences.
#' @param msa_method the method used to create the multiple sequence alignment 
#' @return A list containing:
#' \describe{
#'   \item{Nucleotide_Diversity}{Nucleotide diversity (π).}
#'   \item{Segregating_Sites}{The number of segregating sites.}
#'   \item{Tajima_D}{Tajima's D statistic.}
#' }
#' @details This function converts allele sequences into a DNA alignment format,
#' calculates the specified statistics, and returns them as a list.
#' @examples
#' \dontrun{
#'   stats <- calculate_popgen_stats(c("ATGC", "ATCC", "ATGG"))
#'   print(stats)
#' }
#' @importFrom pegas nuc.div seg.sites tajima.test
#' @export
calculate_popgen_stats <- function(allele_data, msa_method = "Muscle") {

  # Get indices of unique sequences
  unique_seqs <- ! duplicated(allele_data)
  unique_ids <- which(unique_seqs)
  # Extract only unique sequences for alignment
  allele_data_unique <- allele_data[unique_ids]
  # Add names to allow msa() to preserve order
  names(allele_data_unique) <- 1:length(allele_data_unique)
  # Create a mapping from original to unique
  orig2unique_mapping <- match(allele_data, allele_data_unique)

  # Return early if all sequences identical
  if (length(allele_data_unique) == 1) {
    return(list(
      Nucleotide_Diversity = 0, 
      Segregating_Sites = 0, 
      Tajima_D = 0
    ))
  }
  
  # Align unique sequences
  aligned_unique <- msa::msa(
      allele_data_unique, 
      method = msa_method, 
      type = "dna", 
      order = "input"
    ) %>%
    msa::msaConvert(type = "ape::DNAbin")
  # Check order
  if (! identical(names(allele_data_unique), labels(aligned_unique))) {
    stop("Order does not match between alignment input and output")
  }
  # Index aligned sequences with mapping to restore duplicates
  aligned_all <- aligned_unique[orig2unique_mapping, ]
  
  # Compute pop. gen. stats
  nucleotide_diversity <- nuc.div(aligned_all)
  segregating_sites <- length(seg.sites(aligned_all))
  tajima_test <- tajima.test(aligned_all)
 
  # Return the results as a list
  return(list(
    Nucleotide_Diversity = nucleotide_diversity,
    Segregating_Sites = segregating_sites,
    Tajima_D = tajima_test$D, 
    Tajima_D_pval_normal = tajima_test$Pval.normal, 
    Tajima_D_pval_beta = tajima_test$Pval.beta
  ))
}

# Calculate nucleotide diversity, segregating sites, and Tajima's D ------------
#' Calculate Population Genetics Statistics by Target ID
#'
#' Groups allele data by `target_name` and calculates population genetic statistics
#' (nucleotide diversity, number of segregating sites, and Tajima's D) for each group.
#'
#' @param locus_data A data frame containing columns `specimen_name`, `target_name`, and `target_value`.
#' @param msa_method the method used to create the multiple sequence alignment 
#' @return a table of results.
#' Each row corresponds to a `target_name`, and columns include:
#' \describe{
#'   \item{target_name}{The target identifier.}
#'   \item{Nucleotide_Diversity}{Nucleotide diversity (π).}
#'   \item{Segregating_Sites}{The number of segregating sites.}
#'   \item{Tajima_D}{Tajima's D statistic.}
#' }
#' @details This function computes population genetic statistics for each unique `target_name` in the input data,
#' and writes the results to a tab-separated file.
#' @examples
#' \dontrun{
#'   calculate_stats_by_target_name(locus_data)
#' }
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr unnest_wider
#' @importFrom readr write_tsv
#' @export
calculate_stats_by_target_name <- function(locus_data, msa_method = "Muscle") {
  results <- locus_data %>%
    dplyr::group_by(target_name) %>%
    dplyr::summarise(
      stats = list(calculate_popgen_stats(target_value, msa_method))
    ) %>%
    tidyr::unnest_wider(stats)
  return(results)
}


if(!(arg$msa_method %in% c('ClustalW', 'ClustalOmega', 'Muscle'))){
  stop(paste0("--msa_method must be 'ClustalW', 'ClustalOmega', or 'Muscle', not ", arg$msa_method))
}

locus_data = create_locus_data(
  arg$allele_table, 
  specimen_name_col = arg$specimen_name_col, 
  target_name_col = arg$target_name_col, 
  target_value_col = arg$target_value_col
)

# Calculate nuc gens
res = calculate_stats_by_target_name(locus_data)
colnames(res) = tolower(colnames(res))
res %>%
  readr::write_tsv(arg$out)
