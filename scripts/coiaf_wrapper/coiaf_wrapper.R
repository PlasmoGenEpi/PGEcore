#!/usr/bin/env Rscript

# =============================================================================
# COIAF Wrapper Script
# =============================================================================
# 
# This script estimates the complexity of infection (COI) using the COIAF package.
# It processes SNP data and optionally population-level minor allele frequencies
# to estimate COI using both frequency and variant-based methods.
#
# Usage:
#   Rscript coiaf_wrapper.R --snp_data <path> --output <path> [options]
#
# Required arguments:
#   --snp_data    Path to SNP data file (TSV format)
#   --output      Path for output file (TSV format)
#
# Optional arguments:
#   --plmaf       Path to population-level minor allele frequency file (TSV format)
#   --seq_error   Sequencing error rate (default: 0.01)
#   --max_coi     Maximum COI to consider (default: 25)
#   --help        Show this help message
#
# Input file formats:
#   SNP data (required): TSV file with columns:
#     - specimen_id: The specimen ID
#     - snp_name: The SNP name
#     - read_count: The read count
#     - seq_base: The sequence base (A, C, G, T)
#
#   PLMAF data (optional): TSV file with columns:
#     - snp_name: The SNP name
#     - seq_base: The base (A, C, G, T)
#     - plmaf: The population allele frequency
#
# Output:
#   TSV file with columns:
#     - specimen_id: The specimen ID
#     - coi_freq: COI estimate using frequency method
#     - coi_variant: COI estimate using variant method
#
# Example:
#   Rscript coiaf_wrapper.R --snp_data data/snps.tsv --output results/coi_estimates.tsv
#   Rscript coiaf_wrapper.R --snp_data data/snps.tsv --plmaf data/plmaf.tsv --output results/coi_estimates.tsv --seq_error 0.005 --max_coi 30
#
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(coiaf)
  library(dplyr)
  library(tibble)
  library(readr)
  library(optparse)
})

#' Estimate Complexity of Infection (COI) using COIAF
#'
#' @param snp_data A data frame with the following columns:
#' - specimen_id: The specimen ID
#' - snp_name: The SNP name
#' - read_count: The read count
#' - seq_base: the sequence base (A, C, G, T)
#' @param plmaf An optional data frame with the following columns:
#' - snp_name: The SNP name
#' - seq_base: the base (A, C, G, T)
#' - plmaf: the population allele frequency
#' If not provided, the function will calculate the plmaf using the snp_data
#' @param seq_error Sequencing error rate (default: 0.01)
#' @param max_coi Maximum COI to consider (default: 25)
#' @return A data frame with COI estimates for each specimen
run_coiaf <- function(snp_data, plmaf = NULL, seq_error = 0.01, max_coi = 25) {
  cat("Processing SNP data...\n")

  # Process SNP data to calculate within-sample minor allele frequencies
  processed <- snp_data |>
    dplyr::group_by(specimen_id, snp_name) |>
    dplyr::mutate(coverage = sum(read_count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(wsmaf = read_count / coverage) |>
    dplyr::select(specimen_id, snp_name, seq_base, wsmaf, read_count) |>
    tidyr::complete(
      specimen_id, tidyr::nesting(snp_name, seq_base),
      fill = list(wsmaf = 0, read_count = 0)
    )

  # Calculate population-level minor allele frequencies if not provided
  if (is.null(plmaf)) {
    cat("Calculating population-level minor allele frequencies...\n")
    plmaf <- snp_data |>
      dplyr::group_by(snp_name, seq_base) |>
      dplyr::summarize(read_count = sum(read_count), .groups = "drop") |>
      dplyr::group_by(snp_name) |>
      dplyr::mutate(plmaf = read_count / sum(read_count)) |>
      dplyr::arrange(plmaf, .by_group = TRUE) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(snp_name, seq_base, plmaf) |>
      dplyr::filter(plmaf < 1)
  }

  # Filter processed data to only include targets with PLMAF data
  filtered_processed <- processed |>
    dplyr::right_join(plmaf, by = c("snp_name", "seq_base"))

  cat("Estimating COI for", length(unique(filtered_processed$specimen_id)), "specimens...\n")
  
  # Estimate COI using both frequency and variant methods
  coi <- filtered_processed |>
    dplyr::group_by(specimen_id) |>
    dplyr::summarize(
      coi_freq = coiaf::optimize_coi(
        tibble::tibble(wsmaf, plmaf, coverage = read_count),
        data_type = "real", coi_method = "frequency",
        seq_error = seq_error, max_coi = max_coi
      ),
      coi_variant = coiaf::optimize_coi(
        tibble::tibble(wsmaf, plmaf, coverage = read_count),
        data_type = "real", coi_method = "variant",
        seq_error = seq_error, max_coi = max_coi
      ),
      .groups = "drop"
    )
  
  return(coi)
}

#' Validate input data format
#' 
#' @param data Data frame to validate
#' @param required_cols Vector of required column names
#' @param data_name Name of the data for error messages
#' @return TRUE if valid, stops execution if invalid
validate_data <- function(data, required_cols, data_name) {
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in", data_name, ":", paste(missing_cols, collapse = ", ")))
  }
  
  if (nrow(data) == 0) {
    stop(paste(data_name, "is empty"))
  }
  
  return(TRUE)
}

#' Main function to run the COIAF analysis
main <- function() {
  # Define command line options
  option_list <- list(
    make_option(c("--snp_data"), type = "character", default = NULL,
                help = "Path to SNP data file (TSV format) [required]",
                metavar = "FILE"),
    make_option(c("--plmaf"), type = "character", default = NULL,
                help = "Path to population-level minor allele frequency file (TSV format) [optional]",
                metavar = "FILE"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "Path for output file (TSV format) [required]",
                metavar = "FILE"),
    make_option(c("--seq_error"), type = "numeric", default = 0.01,
                help = "Sequencing error rate [default: %default]",
                metavar = "NUM"),
    make_option(c("--max_coi"), type = "integer", default = 25,
                help = "Maximum COI to consider [default: %default]",
                metavar = "NUM")
  )
  
  # Create option parser
  parser <- OptionParser(
    usage = "Usage: %prog --snp_data FILE --output FILE [options]",
    option_list = option_list,
    description = "COIAF Wrapper Script - Estimates complexity of infection (COI) using the COIAF package.",
    epilogue = "Example:\n  %prog --snp_data data/snps.tsv --output results/coi_estimates.tsv\n  %prog --snp_data data/snps.tsv --plmaf data/plmaf.tsv --output results/coi_estimates.tsv --seq_error 0.005 --max_coi 30"
  )
  
  # Parse arguments
  opt <- parse_args(parser, positional_arguments = FALSE)
  
  # Validate required arguments
  if (is.null(opt$snp_data)) {
    cat("Error: --snp_data is required\n")
    print_help(parser)
    quit(status = 1)
  }
  
  if (is.null(opt$output)) {
    cat("Error: --output is required\n")
    print_help(parser)
    quit(status = 1)
  }
  
  # Validate file paths
  if (!file.exists(opt$snp_data)) {
    stop(paste("SNP data file not found:", opt$snp_data))
  }
  
  if (!is.null(opt$plmaf) && !file.exists(opt$plmaf)) {
    stop(paste("PLMAF file not found:", opt$plmaf))
  }
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(opt$output)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Validate parameters
  if (opt$seq_error < 0 || opt$seq_error > 1) {
    stop("seq_error must be between 0 and 1")
  }
  
  if (opt$max_coi < 1) {
    stop("max_coi must be at least 1")
  }
  
  cat("COIAF Wrapper Script\n")
  cat("===================\n")
  cat("SNP data file:", opt$snp_data, "\n")
  if (!is.null(opt$plmaf)) {
    cat("PLMAF file:", opt$plmaf, "\n")
  }
  cat("Output file:", opt$output, "\n")
  cat("Sequencing error rate:", opt$seq_error, "\n")
  cat("Maximum COI:", opt$max_coi, "\n")
  cat("\n")
  
  tryCatch({
    # Read SNP data
    cat("Reading SNP data...\n")
    snp_data <- readr::read_tsv(opt$snp_data, show_col_types = FALSE)
    validate_data(snp_data, c("specimen_id", "snp_name", "read_count", "seq_base"), "SNP data")
    
    # Read PLMAF data if provided
    plmaf_data <- NULL
    if (!is.null(opt$plmaf)) {
      cat("Reading PLMAF data...\n")
      plmaf_data <- readr::read_tsv(opt$plmaf, show_col_types = FALSE)
      validate_data(plmaf_data, c("snp_name", "seq_base", "plmaf"), "PLMAF data")
    }
    
    # Run COIAF analysis
    results <- run_coiaf(
      snp_data = snp_data,
      plmaf = plmaf_data,
      seq_error = opt$seq_error,
      max_coi = opt$max_coi
    )
    
    # Write results
    cat("Writing results to", opt$output, "...\n")
    readr::write_tsv(results, opt$output)
    
    cat("Analysis completed successfully!\n")
    cat("Results saved to:", opt$output, "\n")
    cat("Number of specimens processed:", nrow(results), "\n")
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}

# Run the main function if this script is executed directly
if (!interactive()) {
  main()
}


