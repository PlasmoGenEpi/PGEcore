library(dplyr)
library(tidyr)
library(xlsx)
library(optparse)
library(stringr)
source("SNPModel.R")

# Parse arguments ------------------------------------------------------
opts = list(
  make_option(
    "--aa_calls",
    type = "character", 
    help = str_c(
      "TSV containing aminoacid calls per specimen, with the columns: specimen_id,
      target_id, specimen_id, gene_id, aa_position, ref_aa, and aa."
    )
  )
)

arg <- parse_args(OptionParser(option_list = opts))

# MultiLociBiallelicModel_wrapper functions ----------------------------

#' Create Input for Multi-Loci Biallelic Model
#'
#' This function reads an allele table from a file, processes it to create an input dataset for a multi-loci biallelic model, and validates the resulting data structure.
#'
#' @param input_path A string specifying the file path to the input data. The file should be a tab-separated CSV file with columns including `specimen_id`, `gene_id`, `aa_position`, `ref_aa`, and `aa`.
#'
#' @return A data frame that has been processed to include:
#' \describe{
#'   \item{specimen_id}{The specimen identifiers from the input file.}
#'   \item{identifier}{A unique identifier for each position, combining `gene_id` and `aa_position`.}
#'   \item{value}{An integer value for each entry:
#'     \itemize{
#'       \item \code{0}: All entries for the position have `ref_aa == aa`.
#'       \item \code{1}: All entries for the position have `ref_aa != aa`.
#'       \item \code{2}: Mixed entries, where some have `ref_aa == aa` and others have `ref_aa != aa`.
#'       \item \code{NA}: Positions with more than two distinct amino acids (`aa`).}}
#' }
#'
#' @details
#' The function applies validation rules to ensure that all values in the resulting data frame (except `specimen_id`) are integers in the range \code{0} to \code{2}, or \code{NA}. If any validation checks fail, the function raises an error and stops execution.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' input_path <- "path/to/allele_table.csv"
#' model_input <- create_MultiLociBiallelicModel_input(input_path)
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import validate
#' @export
create_MultiLociBiallelicModel_input <- function(input_path) {
  # Read the allele table
  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")

  MultiLociBiallelicModel_data <- input_data %>%
    mutate(identifier = paste(gene_id, aa_position, sep = "_")) %>%
    group_by(specimen_id, identifier) %>%
    mutate(value = case_when(
      n_distinct(aa) > 2 ~ NA_real_,                      # More than two distinct amino acids
      all(ref_aa == aa) ~ 0,                              # All entries have ref_aa == aa
      all(ref_aa != aa) ~ 1,                              # All entries have ref_aa != aa
      any(ref_aa == aa) & any(ref_aa != aa) ~ 2           # Mixed entries
    )) %>%
    ungroup() %>%
    select(specimen_id, identifier, value) %>%
    distinct() %>%
    pivot_wider(names_from = identifier, values_from = value, values_fill = NA)

  # Generate validation rules to check if all columns (except `specimen_id`) are numeric
  # Validate package couldn't handle a tibble with variable number of columns.
  all_doubles <- all(sapply(MultiLociBiallelicModel_data[-1], is.double))  
  
  # Raise an error if any validations fail
  if (!all_doubles) {
    stop(
      "Analysis object failed one or more validation checks:\n",
      paste("Amino acid table misformatted. Verify your inputs.", collapse = "\n"),
      call. = FALSE
    )
  }
  return(MultiLociBiallelicModel_data)
}

#' Run Multi-Loci Biallelic Model
#'
#' This function runs the multi-loci biallelic model on the provided input data, estimating parameters and allele frequencies, and returns the results along with runtime statistics.
#'
#' @param inputToMLE A data structure prepared as input for the multi-loci biallelic model. This should typically be created using the \code{create_MultiLociBiallelicModel_input} function or similar preprocessing.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{lambda}{The estimated parameter `lambda` from the model.}
#'   \item{plsf_table}{A tibble containing the allele frequencies, with columns:
#'     \itemize{
#'       \item \code{sequence}: The names of sequences or loci.
#'       \item \code{MLBM_frequency}: The estimated allele frequencies for each sequence.}}
#'   \item{runtime}{A named vector of runtime statistics, as returned by \code{system.time}, including user, system, and elapsed time.}
#' }
#'
#' @details
#' The function uses the \code{mle} function to estimate parameters for the model. If the input data has missing values, additional preprocessing to calculate allele frequencies and impute missing data is currently marked as a TODO.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' input_data <- create_MultiLociBiallelicModel_input("path/to/data.csv")
#' results <- run_MultiLociBiallelicModel(input_data)
#' print(results$lambda)
#' print(results$plsf_table)
#' print(results$runtime)
#' }
#'
#' @importFrom tibble tibble
#' @importFrom stats mle
#' @export
run_MultiLociBiallelicModel <- function(inputToMLE) {
  # TODO: If input has missing data calculate allele frequencies and fill it in
  runtime <- system.time({
    est <- mle(inputToMLE, id = TRUE)
  })
  plsf <- est$p
  
  plsf_table <- tibble(
    sequence = dimnames(plsf)[[2]],
    MLBM_frequency = c(plsf)
  )
  return(
    list(
      lambda = est$lambda,
      plsf_table = plsf_table,
      runtime = runtime
    )
  )
}

# Main -----------------------------------------------------------------

# Read the aminoacid calls

#arg$aa_calls = "example_amino_acid_calls.tsv"

MLBM_data <- create_MultiLociBiallelicModel_input(arg$aa_calls)

# Run the Multi-Loci Biallelic Model

MLBM_res <- run_MultiLociBiallelicModel(MLBM_data)

saveRDS(MLBM_res, "moire_res.rds")

