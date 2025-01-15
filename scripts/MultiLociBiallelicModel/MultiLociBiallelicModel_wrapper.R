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
  
  MLBM_data <- input_data %>%
    mutate(identifier = paste(target_id, aa_position, sep = ":")) %>%
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
  all_doubles <- all(sapply(MLBM_data[-1], is.double))  
  
  # Raise an error if any validations fail
  if (!all_doubles) {
    stop(
      "Analysis object failed one or more validation checks:\n",
      paste("Amino acid table misformatted. Verify your inputs.", collapse = "\n"),
      call. = FALSE
    )
  }
  
  staves_data <- input_data %>%
    mutate(
      staves = paste(target_id, aa_position, aa, sep = ":"),
      state = if_else(
        substr(staves, nchar(staves), nchar(staves)) == ref_aa, 
        0, 
        1
      )
    ) %>%
    distinct(staves, .keep_all = TRUE) %>%
    select(staves, state) %>%
    mutate(prestave = sub(":([^:]+)$", "", staves))
  
  
  MLBM_object = list(MLBM_data = MLBM_data, 
                     staves_data = staves_data)
  return(MLBM_object)
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

#' Summarise Multi-Loci Biallelic Model Results
#'
#' This function processes and summarises results from the Multi-Loci Biallelic Model (MLBM). It updates the sequence matrix with corresponding values from the `staves_data`, modifies the sequence representation, and writes the summarised results to a TSV file.
#'
#' @param MLBM_res A list containing the results of the MLBM analysis, including:
#' \describe{
#'   \item{\code{plsf_table}}{A data frame with columns \code{sequence} (encoded as a string of 0s and 1s) and \code{MLBM_frequency} (the frequency of the sequence).}
#' }
#' @param MLBM_object A list containing the MLBM data and additional information, including:
#' \describe{
#'   \item{\code{MLBM_data}}{A data frame where each column corresponds to a specific genetic locus and each row to a specimen.}
#'   \item{\code{staves_data}}{A data frame with columns:
#'     \itemize{
#'       \item \code{staves}: The encoded genetic representation.
#'       \item \code{state}: The state (0 or 1) of the representation.
#'       \item \code{prestave}: The genetic locus identifier without the encoded suffix.}}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Converts the \code{sequence} column in \code{MLBM_res$plsf_table} into a matrix, splitting each sequence into individual loci.
#'   \item Replaces 0s and 1s in the sequence matrix with corresponding genetic representations from \code{MLBM_object$staves_data}.
#'   \item Updates the \code{sequence} column in \code{MLBM_res$plsf_table} by concatenating all loci values row-wise, separated by semicolons (\code{;}).
#'   \item Renames columns in \code{MLBM_res$plsf_table} to \code{seq} and \code{freq} for clarity.
#'   \item Writes the summarised table to a TSV file named \code{"MLBM_summary.tsv"} in the working directory.
#' }
#'
#' @return None. The function modifies the input \code{MLBM_res} and writes the summarised results to a file.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' summarise_MLBM_results(MLBM_res, MLBM_object)
#' }
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#' @export
summarise_MLBM_results <- function(MLBM_res, MLBM_object) {
  sequence_matrix <- MLBM_res$plsf_table %>%
    mutate(sequence_split = strsplit(sequence, "")) %>%  # Split each sequence into characters
    unnest_wider(sequence_split, names_sep = "_") %>%   # Convert the list of characters to columns
    select(-MLBM_frequency) %>%                         # Remove the frequency column
    rename_with(~ colnames(MLBM_object$MLBM_data)[-1], starts_with("sequence_split")) %>% # Rename columns
    as.matrix() 
  
  # Convert sequence_matrix to a data frame for easier manipulation
  sequence_df <- as.data.frame(sequence_matrix, stringsAsFactors = FALSE)
  columns_to_replace <- colnames(sequence_df)[-1]  # Exclude the "sequence" column
  
  # Replace 0s and 1s in sequence_matrix[,-1] based on staves_data
  for (col in columns_to_replace) {
    matching_rows <- MLBM_object$staves_data %>%
      filter(prestave == col)  # Match column name with prestave
    for(i in 1:length(sequence_df[[col]])) {
      sequence_df[[col]][i] = matching_rows$staves[matching_rows$state == sequence_df[[col]][i]]
    }
  }
  
  # Convert the modified data frame back to a matrix if needed
  sequence_matrix_updated <- as.matrix(sequence_df)
  
  MLBM_res$plsf_table$sequence <- apply(sequence_df[-1], 1, paste, collapse = ";")
  
  MLBM_res$plsf_table <- MLBM_res$plsf_table %>%
    rename(seq = sequence, 
           freq = MLBM_frequency)
  
  readr::write_tsv(MLBM_res$plsf_table, "MLBM_summary.tsv")
}

# Main -----------------------------------------------------------------

# Read the aminoacid calls
#arg$aa_calls = "example_amino_acid_calls.tsv"

MLBM_object <- create_MultiLociBiallelicModel_input(arg$aa_calls)

# Run the Multi-Loci Biallelic Model
MLBM_res <- run_MultiLociBiallelicModel(MLBM_object$MLBM_data)

# Generate summary
summarise_MLBM_results(MLBM_res, MLBM_object)

