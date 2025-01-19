library(optparse)
library(stringr)
library(dplyr)
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
  ),
  make_option(
    "--loci_group_table",
    type = "character", 
    help = str_c(
      "TSV containing loci groups, with the columns: group_id, gene_id, 
      aa_position."
    )
  )
)

arg <- parse_args(OptionParser(option_list = opts))

# MultiLociBiallelicModel_wrapper functions ----------------------------

#' Create Multi-Loci Biallelic Model Input
#'
#' This function processes allele data and generates the input required for multi-loci biallelic analysis.
#' It validates input data, processes amino acid information, matches loci groups, and creates grouped
#' subtables for downstream analysis.
#'
#' @param input_path A string specifying the path to the input allele table (in tab-delimited format).
#' @param loci_group A string specifying the path to the loci group file (in tab-delimited format).
#'
#' @return A list (`MLBM_object`) containing:
#' \describe{
#'   \item{`MLBM_data`}{A processed tibble containing the input allele data in a wide format.}
#'   \item{`staves_data`}{A tibble mapping each allele to its state (reference or alternative).}
#'   \item{`by_group_table`}{A list of tibbles, where each tibble corresponds to a group of loci and contains the associated allele data.}
#'   \item{`groups`}{A character vector of unique `group_id` values.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Reads the allele table from `input_path` and processes it to create a wide-format tibble.
#'   \item Validates the allele data to ensure all loci columns (excluding `specimen_id`) are numeric.
#'   \item Processes the loci group file from `loci_group` and validates it against predefined rules.
#'   \item Matches loci from the allele table to their corresponding groups, creating subtables for each group.
#'   \item Generates a `MLBM_object` containing all processed data and grouped subtables.
#' }
#'
#' @examples
#' # Define file paths
#' input_path <- "allele_table.tsv"
#' loci_group <- "loci_group.tsv"
#'
#' # Generate the Multi-Loci Biallelic Model input
#' MLBM_object <- create_MultiLociBiallelicModel_input(input_path, loci_group)
#'
#' @import dplyr
#' @import tidyr
#' @import validate
#' @export
create_MultiLociBiallelicModel_input <- function(input_path, loci_group) {
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
    tidyr::pivot_wider(names_from = identifier, values_from = value, values_fill = NA)
  
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
  
  loci_groups = read.csv(loci_group, sep = "\t")
  
  rules <- validate::validator(
    # Data columns
    is.character(group_id),
    is.character(gene_id),
    is.integer(aa_position),
    
    # Non-missing values
    !is.na(group_id),
    !is.na(gene_id),
    !is.na(aa_position)
  )
  
  # Confront the analysis_object with validation rules
  print("Confronting input data with validation rules")
  fails <- validate::confront(loci_groups, rules, raise = "all") %>%
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
  
  #Matching groups to create lists of loci
  match_group = input_data %>% 
    select(target_id, gene_id, aa_position) %>%
    mutate(identifier = paste(target_id, aa_position, sep = ":")) 
  
  merged_data <- match_group %>%
    inner_join(loci_groups, by = c("gene_id", "aa_position"), relationship = "many-to-many")
  
  grouped_list <- merged_data %>%
    group_by(group_id) %>%
    summarise(identifiers = list(unique(identifier)), .groups = "drop")
  
  result_list <- setNames(grouped_list$identifiers, grouped_list$group_id)
  
  # Generate subtables
  result_tables <- lapply(names(result_list), function(group_name) {
    columns <- c("specimen_id", result_list[[group_name]])
    subtable <- MLBM_data %>% select(all_of(columns))
    return(subtable)
  })
  
  # Name the list elements
  names(result_tables) <- names(result_list)
  MLBM_object$by_group_table <- result_tables
  MLBM_object$groups = unique(loci_groups$group_id)
  
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
#' This function processes and summarises results from a Multi-Loci Biallelic Model (MLBM) analysis.
#' It modifies the sequence representation in `MLBM_res` based on `MLBM_object` and the provided `group_name`,
#' updating the `plsf_table` with new sequence variants and frequencies.
#'
#' @param MLBM_res A list containing the results of the MLBM analysis, including the `plsf_table`.
#' @param MLBM_object A list containing the input data and metadata used in the MLBM analysis, including `by_group_table` and `staves_data`.
#' @param group_name A string specifying the name of the group to process, corresponding to a key in `MLBM_object$by_group_table`.
#'
#' @return The modified `MLBM_res` object with an updated `plsf_table` that includes:
#' \describe{
#'   \item{`group_id`}{A column indicating the group name.}
#'   \item{`variant`}{A column representing the reconstructed sequence variants.}
#'   \item{`freq`}{A column indicating the frequency of each variant.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Splits sequences in `MLBM_res$plsf_table` into individual characters and converts them to a wide format.
#'   \item Matches and replaces sequence states using `staves_data` in `MLBM_object`.
#'   \item Reconstructs updated sequences by combining individual character columns into a single string.
#'   \item Updates the `plsf_table` with new `group_id`, `variant`, and `freq` columns.
#' }
#'
#' @examples
#' # Example usage
#' summarised_results <- summarise_MLBM_results(MLBM_res, MLBM_object, "pfdhfr_pfdhps")
#'
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @export
summarise_MLBM_results <- function(MLBM_res, MLBM_object, group_name) {
  sequence_matrix <- MLBM_res$plsf_table %>%
    mutate(sequence_split = strsplit(sequence, "")) %>%  # Split each sequence into characters
    tidyr::unnest_wider(sequence_split, names_sep = "_") %>%   # Convert the list of characters to columns
    select(-MLBM_frequency) %>%                         # Remove the frequency column
    rename_with(~ colnames(MLBM_object$by_group_table[[group_name]])[-1], starts_with("sequence_split")) %>% # Rename columns
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
  
  # Add group_id column and rename columns
  MLBM_res$plsf_table <- MLBM_res$plsf_table %>%
    mutate(group_id = group_name, .before = 1) %>% # Add group_id at the beginning
    rename(variant = sequence, freq = MLBM_frequency) # Rename columns
}

# Main -----------------------------------------------------------------

# Read the aminoacid calls

MLBM_object <- create_MultiLociBiallelicModel_input(arg$aa_calls, arg$loci_group)

MLBM_res = list()
# Iterate over each table by groups of loci
for (group_name in MLBM_object$groups) {
  print(paste0("Running MultiLociBiallelicModel on ", group_name))
  group_table <- MLBM_object$by_group_table[[group_name]]
  MLBM_tmp <- run_MultiLociBiallelicModel(group_table)
  MLBM_res[[group_name]] <- summarise_MLBM_results(MLBM_tmp, MLBM_object, group_name)
}

readr::write_tsv(bind_rows(MLBM_res), "MLBM_summary.tsv")