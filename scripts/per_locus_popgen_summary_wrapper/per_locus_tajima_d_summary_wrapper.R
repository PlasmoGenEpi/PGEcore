#!/usr/bin/env Rscript

library("optparse")
library("stringr")
library("pegas")

# Parse arguments ------------------------------------------------------
opts = list(
  make_option(
    "--allele_table",
    help = str_c(
      "TSV containing allele present/absent per specimen, with the
       columns: specimen_id, target_id, seq"
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




# locus counting functions -----------------------------------------------------

# create locus data ------------------------------------------------------------
#' Create Locus Data
#'
#' Reads an input file, validates its format, and creates a locus data frame
#' containing sample IDs, target IDs, and allele sequences.
#'
#' @param input_path A string specifying the path to the input file. The file should be tab-separated
#' and contain columns for `specimen_id`, `target_id`, and `seq`.
#' @return A data frame with columns `sample_id`, `target_id`, and `allele`.
#' @details This function reads the input data from a file, validates the format using
#' predefined rules (ensuring all values are non-missing and of the correct type), and
#' returns a cleaned data frame suitable for further analysis.
#' @examples
#' \dontrun{
#'   locus_data <- create_locus_data("path/to/input_file.tsv")
#' }
#' @importFrom dplyr select rename filter
#' @importFrom validate validator confront summary
#' @importFrom stringr str_c
#' @export
create_locus_data <- function(input_path) {

  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")
  locus_data <- input_data |>
    dplyr::select(specimen_id, target_id, seq) |> 
    dplyr::rename(sample_id = specimen_id, allele = seq)
  
  print("Validating input format")
  rules <- validate::validator(
    # Data columns
    is.character(sample_id),
    is.character(target_id),
    is.character(allele),
  
    # Non-missing values
    !is.na(sample_id),
    !is.na(target_id),
    !is.na(allele)
  )
  
  # Confront the analysis_object with validation rules
  print("Confronting input data with validation rules")
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
  
  print("Returning Locus data")
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
  alignment <- msa::msa(allele_data, method = msa_method, type = "dna")
  aligned_sequences <- msa::msaConvert(alignment, type = "ape::DNAbin")
  
  nucleotide_diversity <- nuc.div(aligned_sequences)
  segregating_sites <- length(seg.sites(aligned_sequences))
  if(segregating_sites == 0){
    tajima_test = list(D = 0)
  } else {
    tajima_test <- tajima.test(aligned_sequences)
  }
 
  # Return the results as a list
  return(list(
    Nucleotide_Diversity = nucleotide_diversity,
    Segregating_Sites = segregating_sites,
    Tajima_D = tajima_test$D
  ))
}

# Calculate nucleotide diversity, segregating sites, and Tajima's D ------------
#' Calculate Population Genetics Statistics by Target ID
#'
#' Groups allele data by `target_id` and calculates population genetic statistics
#' (nucleotide diversity, number of segregating sites, and Tajima's D) for each group.
#'
#' @param locus_data A data frame containing columns `sample_id`, `target_id`, and `allele`.
#' @param msa_method the method used to create the multiple sequence alignment 
#' @return a table of results.
#' Each row corresponds to a `target_id`, and columns include:
#' \describe{
#'   \item{target_id}{The target identifier.}
#'   \item{Nucleotide_Diversity}{Nucleotide diversity (π).}
#'   \item{Segregating_Sites}{The number of segregating sites.}
#'   \item{Tajima_D}{Tajima's D statistic.}
#' }
#' @details This function computes population genetic statistics for each unique `target_id` in the input data,
#' and writes the results to a tab-separated file.
#' @examples
#' \dontrun{
#'   calculate_stats_by_target_id(locus_data)
#' }
#' @importFrom dplyr group_by summarise
#' @importFrom tidyr unnest_wider
#' @importFrom readr write_tsv
#' @export
calculate_stats_by_target_id <- function(locus_data, msa_method = "Muscle") {
  results <- locus_data %>%
    dplyr::group_by(target_id) %>%
    dplyr::summarise(
      stats = list(calculate_popgen_stats(allele, msa_method))
    ) %>%
    tidyr::unnest_wider(stats)
  return(results)
}

# Main function ------------------------------------------------------
# Load allele table/locus data
#arg$allele_table = "/Users/jar4142/Desktop/PGEcore/data/example2_allele_table.tsv"
arg <- parse_args(OptionParser(option_list = opts))

if(!(arg$msa_method %in% c('ClustalW', 'ClustalOmega', 'Muscle'))){
  stop(paste0("--msa_method must be 'ClustalW', 'ClustalOmega', or 'Muscle', not ", arg$msa_method))
}

locus_data = create_locus_data(arg$allele_table)

# Calculate nuc gens
res = calculate_stats_by_target_id(locus_data)
colnames(res) = tolower(colnames(res))
readr::write_tsv(res, arg$out)

