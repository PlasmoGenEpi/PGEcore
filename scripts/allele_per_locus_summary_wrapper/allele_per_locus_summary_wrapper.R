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

create_locus_data <- function(input_path) {

  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")
  locus_data <- input_data |>
    dplyr::select(specimen_id, target_id, seq) |> 
    dplyr::rename(sample_id = specimen_id, locus = target_id, allele = seq)
  
  print("Validating input format")
  rules <- validate::validator(
    # Data columns
    is.character(sample_id),
    is.character(locus),
    is.character(allele),
  
    # Non-missing values
    !is.na(sample_id),
    !is.na(locus),
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

# Function to calculate allele metrics for each target_id
summarize_allele_table <- function(locus_data) {
  # Split sequences into individual alleles and calculate metrics per target_id
  allele_summary_by_target <- locus_data %>%
    group_by(locus) %>%
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
#arg$allele_table <- "/Users/jar4142/Desktop/PGEcore/data/example2_allele_table.tsv"

# Load allele table/locus data
locus_data = create_locus_data(arg$allele_table)

# Generate locus summary
summarize_allele_table(locus_data)

