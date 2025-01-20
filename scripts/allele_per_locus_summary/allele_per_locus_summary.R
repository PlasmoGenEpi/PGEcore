library(dplyr)
library(tidyr)
library(optparse)
library(stringr)

# Definitions:
# By target_id
# Total Allele Count: The total number of alleles per target_id (length(alleles)).
# Unique Allele Count: The count of unique alleles (length(unique(alleles))).
#	Allele Singlets: The number of alleles that appear only once sum(table(seq) == 1)

# Parse arguments ------------------------------------------------------
opts = list(
  make_option(
    "--allele_table",
    help = str_c(
      "TSV containing allele present/absent per specimen, with the
       columns: specimen_id, target_id, seq"
    )
  )
)

arg <- parse_args(OptionParser(option_list = opts))

# locus counting functions ----------------------------------------------------

#' Create Locus Data from Input File
#'
#' This function reads input data from a specified file, validates its format, and extracts the required locus information.
#'
#' @param input_path A string specifying the file path to the input data. The file should be a tab-separated CSV containing the columns:
#' \describe{
#'   \item{\code{specimen_id}}{Identifier for each specimen.}
#'   \item{\code{target_id}}{Identifier for each genetic target.}
#'   \item{\code{seq}}{Sequence data representing alleles.}
#' }
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{\code{sample_id}}{Renamed from \code{specimen_id}, representing the specimen identifier.}
#'   \item{\code{target_id}}{The genetic target identifier.}
#'   \item{\code{allele}}{Renamed from \code{seq}, representing allele sequence data.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Reads the input data file and selects the relevant columns (\code{specimen_id}, \code{target_id}, and \code{seq}).
#'   \item Renames the columns to \code{sample_id}, \code{target_id}, and \code{allele} for consistency.
#'   \item Validates the data format to ensure all columns are of type character and contain no missing values.
#'   \item Raises an error if any of the validation checks fail, providing details on the failed expressions.
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' input_path <- "example_locus_data.tsv"
#' locus_data <- create_locus_data(input_path)
#' print(locus_data)
#' }
#'
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

# Calculate allele metrics for each target_id --------------------------------
#' Summarize Allele Metrics by Target ID
#'
#' This function calculates summary metrics for alleles grouped by `target_id` and writes the results to a TSV file.
#'
#' @param locus_data A data frame containing locus data with at least the following columns:
#' \describe{
#'   \item{\code{target_id}}{The genetic target identifier.}
#'   \item{\code{allele}}{The sequence data representing alleles.}
#' }
#'
#' @return None. The function writes a summary table to a file named \code{"allele_summary_by_target.tsv"} in the working directory.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Groups the data by \code{target_id}.
#'   \item Calculates the following metrics for each \code{target_id}:
#'     \itemize{
#'       \item \code{total_allele_count}: The total number of alleles (\code{length(allele)}).
#'       \item \code{unique_allele_count}: The number of unique alleles (\code{length(unique(allele))}).
#'       \item \code{allele_singlets}: The count of alleles that appear only once (\code{sum(table(allele) == 1)}).
#'     }
#'   \item Writes the summary table to a TSV file.
#' }
#'
#' The output table includes the following columns:
#' \describe{
#'   \item{\code{target_id}}{The genetic target identifier.}
#'   \item{\code{total_allele_count}}{The total number of alleles for the target.}
#'   \item{\code{unique_allele_count}}{The number of unique alleles for the target.}
#'   \item{\code{allele_singlets}}{The number of alleles that appear only once for the target.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' allele_summary <- summarize_allele_table(locus_data)
#' }
#'
#' @importFrom dplyr group_by summarize ungroup
#' @importFrom readr write_tsv
#' @export

summarize_allele_table <- function(locus_data) {
  # Split sequences into individual alleles and calculate metrics per target_id
  allele_summary_by_target <- locus_data %>%
    group_by(target_id) %>%
    summarize(
      total_allele_count = length(allele),                    # Total allele count
      unique_allele_count = length(unique(allele)), # Unique allele count
      allele_singlets = sum(table(allele) == 1)    # Number of singlets
    ) %>%
    ungroup()
  
  # Write summary to files
  readr::write_tsv(allele_summary_by_target, "allele_summary_by_target.tsv")
}

# Main-----------------------------------------------------------------

# Load allele table/locus data
locus_data = create_locus_data(arg$allele_table)

# Generate locus summary
summarize_allele_table(locus_data)

