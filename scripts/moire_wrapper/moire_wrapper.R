library(tibble)
library(dplyr)
library(moire)
library(tidyr)
library(optparse)
library(stringr)
library(validate)
library(checkmate)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table",
    help = str_c(
       "TSV containing allele present/absent per specimen, with the
       columns: specimen_id, target_id, seq"
    )
  ),
  make_option(
    "--allow_relatedness",
    type = "logical",
    default = TRUE,
    help = str_c(
      "Logical flag to allow relatedness within samples.",
      "Set to TRUE by default"
    )
  ),
  make_option(
    "--burnin",
    type = "integer",
    default = 10000,
    help = str_c(
      "Number of burn-in iterations for the MCMC sampler.",
      "Default is set to 10000 iterations."
    )
  ),
  make_option(
    "--samples_per_chain",
    type = "integer",
    default = 1000,
    help = str_c(
      "Number of samples per chain for the MCMC sampler.",
      "Default is set to 1000 samples."
    )
  ),
  make_option(
    "--verbose",
    type = "logical",
    default = FALSE,
    help = str_c(
      "Logical flag to enable verbose output during processing.",
      "Set to FALSE by default."
    )
  ),
  make_option(
    "--eps_pos_alpha",
    type = "integer",
    default = 1,
    help = str_c(
      "Alpha parameter for the positive error rate prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--eps_pos_beta",
    type = "integer",
    default = 1,
    help = str_c(
      "Beta parameter for the positive error rate prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--eps_neg_alpha",
    type = "integer",
    default = 1,
    help = str_c(
      "Alpha parameter for the negative error rate prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--eps_neg_beta",
    type = "integer",
    default = 1,
    help = str_c(
      "Beta parameter for the negative error rate prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--r_alpha",
    type = "integer",
    default = 1,
    help = str_c(
      "Alpha parameter for the relatedness prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--r_beta",
    type = "integer",
    default = 1,
    help = str_c(
      "Beta parameter for the relatedness prior.",
      "Default set to 1."
    )
  ),
  make_option(
    "--mean_coi_shape",
    type = "double",
    default = 0.1,
    help = str_c(
      "Shape parameter for the mean complexity of infection prior.",
      "Default set to 0.1."
    )
  ),
  make_option(
    "--mean_coi_scale",
    type = "integer",
    default = 10,
    help = str_c(
      "Scale parameter for the mean complexity of infection prior.",
      "Default set to 10"
    )
  ),
  make_option(
    "--max_eps_pos",
    type = "integer",
    default = 2,
    help = str_c(
      "Maximum value for the positive error rate.",
      "Default set to 2."
    )
  ),
  make_option(
    "--max_eps_neg",
    type = "integer",
    default = 2,
    help = str_c(
      "Maximum value for the negative error rate.",
      "Default set to 2."
    )
  ),
  make_option(
    "--record_latent_genotypes",
    type = "logical",
    default = FALSE,
    help = str_c(
      "Logical flag to record latent genotypes in the output.",
      "Default set to FALSE."
    )
  ),
  make_option(
    "--pt_chains",
    type = "integer",
    default = 1,
    help = str_c(
      "Number of chains to run in the parallel tempering process.",
      "Default set to 1."
    )
  ),
  make_option(
    "--pt_grad_lower",
    type = "double",
    default = 0,
    help = str_c(
      "Lower bound for the gradient step size in the parallel tempering process.",
      "Default set to 0."
    )
  ),
  make_option(
    "--pt_num_threads",
    type = "integer",
    default = 1,
    help = str_c(
      "Number of threads to use in the parallel tempering process.",
      "Default set to 1."
    )
  ),
  make_option(
    "--adapt_temp",
    type = "logical",
    default = TRUE,
    help = str_c(
      "Logical flag to enable adaptive temperature in the MCMC process.",
      "Default set to TRUE."
    )
  ),
  make_option(
    "--max_runtime",
    type = "double",
    default = Inf,
    help = str_c(
      "Maximum runtime for the MCMC process.",
      "Default set to Inf."
    )
  ),
  make_option(
    "--coi_summary",
    default = "coi_summary.tsv",
    help = str_c(
      "TSV file for COI (complexity of infection) summary. ",
      "This file will contain summarized COI information for each sample."
    )
  ),
  make_option(
    "--he_summary",
    default = "he_summary.tsv",
    help = str_c(
      "TSV file for He (expected heterozygosity) summary. ",
      "This file will provide He values calculated across all loci."
    )
  ),
  make_option(
    "--allele_freq_summary",
    default = "allele_freq_summary.tsv",
    help = str_c(
      "TSV file for allele frequency summary. ",
      "This file will include allele frequencies for all loci."
    )
  ),
  make_option(
    "--relatedness_summary",
    default = "relatedness_summary.tsv",
    help = str_c(
      "TSV file for relatedness summary. ",
      "This file will contain pairwise relatedness estimates for all samples."
    )
  ),
  make_option(
    "--effective_coi_summary",
    default = "effective_coi_summary.tsv",
    help = str_c(
      "TSV file for effective COI summary. ",
      "This file will provide the effective complexity of infection estimates."
    )
  )
)

arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table.tsv"
  arg$coi_summary <- "coi_summary.tsv"
  arg$he_summary <- "he_summary.tsv"
  arg$allele_freq_summary <- "allele_freq_summary.tsv"
  arg$relatedness_summary <- "relatedness_summary.tsv"
  arg$effective_coi_summary <- "effective_coi_summary.tsv"
}

# moire_wrapper functions ------------------------------------------------------

#' Create Moire Input Object
#'
#' This function reads input data, validates its format, and constructs a Moire input object
#' containing data and parameters required for downstream analysis.
#'
#' @param input_path Character. Path to the input data file (CSV format, tab-separated).
#' @param allow_relatedness Logical. Whether to allow relatedness in the analysis.
#' @param burnin Numeric. Number of burn-in iterations for the MCMC algorithm.
#' @param samples_per_chain Numeric. Number of samples per chain in the MCMC algorithm.
#' @param verbose Logical. Whether to display detailed messages during execution.
#' @param eps_pos_alpha Numeric. Alpha parameter for the positive error rate prior.
#' @param eps_pos_beta Numeric. Beta parameter for the positive error rate prior.
#' @param eps_neg_alpha Numeric. Alpha parameter for the negative error rate prior.
#' @param eps_neg_beta Numeric. Beta parameter for the negative error rate prior.
#' @param r_alpha Numeric. Alpha parameter for the relatedness prior.
#' @param r_beta Numeric. Beta parameter for the relatedness prior.
#' @param mean_coi_shape Numeric. Shape parameter for the COI distribution prior.
#' @param mean_coi_scale Numeric. Scale parameter for the COI distribution prior.
#' @param max_eps_pos Numeric. Maximum allowable positive error rate.
#' @param max_eps_neg Numeric. Maximum allowable negative error rate.
#' @param record_latent_genotypes Logical. Whether to record latent genotypes during the analysis.
#' @param pt_chains Numeric. Number of chains for parallel tempering.
#' @param pt_grad_lower Numeric. Lower bound for the gradient step size in the parallel tempering process.
#' @param pt_num_threads Numeric. Number of threads for parallel tempering computations.
#' @param adapt_temp Logical. Whether to adapt the temperature during parallel tempering.
#' @param max_runtime Numeric. Maximum runtime allowed for the MCMC algorithm (in seconds).
#'
#' @return A list containing the Moire data and parameters, ready for downstream analysis.
#' @examples
#' \dontrun{
#' moire_input <- create_moire_input(
#'   input_path = "data.csv", allow_relatedness = TRUE, burnin = 1000,
#'   samples_per_chain = 5000, verbose = TRUE, eps_pos_alpha = 2,
#'   eps_pos_beta = 5, eps_neg_alpha = 2, eps_neg_beta = 5, r_alpha = 1,
#'   r_beta = 1, mean_coi_shape = 2, mean_coi_scale = 0.5, max_eps_pos = 0.1,
#'   max_eps_neg = 0.1, record_latent_genotypes = FALSE, pt_chains = 40,
#'   pt_grad_lower = 0.5, pt_num_threads = 20,
#'   adapt_temp = TRUE, max_runtime = 3600
#' )
#' }
#' @export

create_moire_input <- function(input_path, allow_relatedness, burnin,
                               samples_per_chain, verbose, eps_pos_alpha,
                               eps_pos_beta, eps_neg_alpha, eps_neg_beta,
                               r_alpha, r_beta, mean_coi_shape, mean_coi_scale,
                               max_eps_pos, max_eps_neg, record_latent_genotypes,
                               pt_chains, pt_grad_lower,
                               pt_num_threads, adapt_temp, max_runtime) {
  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")
  moire_data <- input_data |>
    dplyr::select(specimen_id, target_id, seq) |>
    dplyr::rename(sample_id = specimen_id, locus = target_id, allele = seq)

  print("Creating Moire object")
  # Create a list containing selected data and parameters

  if (pt_chains > 1) {
    pt_chains <- seq(
      from = pt_grad_lower, to = 1, length.out = pt_chains
    )
  } else {
    pt_chains <- 1
  }

  moire_object <- list(
    moire_data = moire_data,
    moire_parameters = list(
      allow_relatedness = allow_relatedness,
      burnin = burnin,
      samples_per_chain = samples_per_chain,
      verbose = verbose,
      eps_pos_alpha = eps_pos_alpha,
      eps_pos_beta = eps_pos_beta,
      eps_neg_alpha = eps_neg_alpha,
      eps_neg_beta = eps_neg_beta,
      r_alpha = r_alpha,
      r_beta = r_beta,
      mean_coi_shape = mean_coi_shape,
      mean_coi_scale = mean_coi_scale,
      max_eps_pos = max_eps_pos,
      max_eps_neg = max_eps_neg,
      record_latent_genotypes = record_latent_genotypes,
      pt_chains = pt_chains,
      pt_num_threads = pt_num_threads,
      adapt_temp = adapt_temp,
      max_runtime = max_runtime
    )
  )

  # Validate input format
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

  # Validate parameters
  assert_logical(moire_object$moire_parameters$allow_relatedness, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$burnin, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$samples_per_chain, any.missing = FALSE, len = 1)
  assert_logical(moire_object$moire_parameters$verbose, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$eps_pos_alpha, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$eps_pos_beta, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$eps_neg_alpha, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$eps_neg_beta, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$r_alpha, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$r_beta, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$mean_coi_shape, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$mean_coi_scale, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$max_eps_pos, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$max_eps_neg, any.missing = FALSE, len = 1)
  assert_logical(moire_object$moire_parameters$record_latent_genotypes, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$pt_num_threads, any.missing = FALSE, len = 1)
  assert_logical(moire_object$moire_parameters$adapt_temp, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$max_runtime, any.missing = FALSE, len = 1)

  # Confront the analysis_object with validation rules
  print("Confronting input data with validation rules")
  fails <- validate::confront(moire_object$moire_data, rules, raise = "all") %>%
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

  print("Returning Moire object")
  return(moire_object)
}

#' Run Moire Analysis
#'
#' This function performs MCMC-based analysis using the specified Moire input object.
#'
#' @param moire_object List. The Moire input object created by `create_Moire_input`.
#'
#' @return A list containing the results of the MCMC analysis, including posterior estimates.
#' @examples
#' \dontrun{
#' results <- run_Moire(moire_input)
#' }
#' @export

run_moire <- function(moire_object) {

  moire_data <- moire_object$moire_data |>
    moire::load_long_form_data()

  moire_parameters <- moire_object$moire_parameters

  runtime <- system.time({
    # Run the MCMC function
    moire_res <- moire::run_mcmc(
      moire_data, moire_data$is_missing,
      allow_relatedness = moire_parameters$allow_relatedness,
      burnin = moire_parameters$burnin,
      samples_per_chain = moire_parameters$samples_per_chain,
      verbose = moire_parameters$verbose,
      eps_pos_alpha = moire_parameters$eps_pos_alpha,
      eps_pos_beta = moire_parameters$eps_pos_beta,
      eps_neg_alpha = moire_parameters$eps_neg_alpha,
      eps_neg_beta = moire_parameters$eps_neg_beta,
      r_alpha = moire_parameters$r_alpha,
      r_beta = moire_parameters$r_beta,
      mean_coi_shape = moire_parameters$mean_coi_shape,
      mean_coi_scale = moire_parameters$mean_coi_scale,
      max_eps_pos = moire_parameters$max_eps_pos,
      max_eps_neg = moire_parameters$max_eps_neg,
      record_latent_genotypes = moire_parameters$record_latent_genotypes,
      pt_chains = moire_parameters$pt_chains,
      pt_num_threads = moire_parameters$pt_num_threads,
      adapt_temp = moire_parameters$adapt_temp,
      max_runtime = moire_parameters$max_runtime
    )
  })

  return(moire_res)
}

#' Summarize and Write MoIre Analysis Results
#'
#' This function summarizes key metrics from MCMC results of a MoIre analysis and writes the summaries to specified output files. Missing loci are handled and added to the summaries with default values.
#'
#' @param moire_object A list containing MoIre data, including `moire_data` with columns such as `locus`, `allele`, and `sample_id`.
#' @param mcmc_results The MCMC results from the MoIre analysis used to calculate summaries for COI, He, allele frequencies, relatedness, and effective COI.
#' @param coi_summary_o A string specifying the output file path for the COI summary TSV file.
#' @param he_summary_o A string specifying the output file path for the He summary TSV file.
#' @param allele_freq_summary_o A string specifying the output file path for the allele frequency summary TSV file.
#' @param relatedness_summary_o A string specifying the output file path for the relatedness summary TSV file.
#' @param effective_coi_summary_o A string specifying the output file path for the effective COI summary TSV file.
#'
#' @return None. The function writes the updated summaries to the specified output files.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Summarizes COI (complexity of infection), He (expected heterozygosity), allele frequencies, relatedness, and effective COI from the MCMC results using the appropriate MoIre functions.
#'   \item Identifies loci missing from the allele frequency summary and adds them back with default frequency values.
#'   \item Updates the He summary with the total sample count for each target.
#'   \item Ensures consistency between MoIre data loci and the summaries, validating the presence of all loci.
#'   \item Writes the finalized summaries to the specified file paths in TSV format.
#' }
#'
#' @examples
#' # Example usage:
#' summarize_and_write_results(
#'   moire_object = moire_object,
#'   mcmc_results = mcmc_results,
#'   coi_summary_o = "coi_summary.tsv",
#'   he_summary_o = "he_summary.tsv",
#'   allele_freq_summary_o = "allele_freq_summary.tsv",
#'   relatedness_summary_o = "relatedness_summary.tsv",
#'   effective_coi_summary_o = "effective_coi_summary.tsv"
#' )
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#' @importFrom moire summarize_coi summarize_he summarize_allele_freqs summarize_relatedness summarize_effective_coi
#' @export

summarize_and_write_results <- function(moire_object, mcmc_results, coi_summary_o, he_summary_o, allele_freq_summary_o, relatedness_summary_o, effective_coi_summary_o) {

  # Summarize statistics
  coi_summary <- moire::summarize_coi(mcmc_results) %>%
    rename(specimen_id = sample_id,
           coi = post_coi_mean)
  he_summary <- moire::summarize_he(mcmc_results) %>%
    rename(target_id = locus,
           he = post_stat_mean)
  allele_freq_summary <- moire::summarize_allele_freqs(mcmc_results) %>%
    rename(target_id = locus,
           freq = post_allele_freqs_mean,
           seq = allele)
  relatedness_summary <- moire::summarize_relatedness(mcmc_results) %>%
    rename(specimen_id = sample_id,
           within_host_rel = post_relatedness_mean)
  effective_coi_summary <- moire::summarize_effective_coi(mcmc_results) %>%
    rename(specimen_id = sample_id,
           ecoi = post_effective_coi_mean)

  # Moire removes loci with only 1 allele.
  # Add removed loci to allele_freq_summary. Add freq 1 by default.
  missing_target_ids <- unique(moire_object$moire_data$locus[
    !moire_object$moire_data$locus %in% allele_freq_summary$target_id])
  present_target_ids <- unique(moire_object$moire_data$locus[
    moire_object$moire_data$locus %in% allele_freq_summary$target_id])

  assert(length(unique(moire_object$moire_data$locus)) ==
           length(missing_target_ids) + length(present_target_ids))
  assert(all(sort(present_target_ids) == sort(unique(allele_freq_summary$target_id))))

  one_allele_loci <- moire_object$moire_data[
    moire_object$moire_data$locus %in% missing_target_ids,] %>%
    select(-sample_id) %>%            # Remove the `sample_id` column
    rename(target_id = locus,
           seq = allele) %>%
    distinct() %>%                    # Remove duplicated rows
    mutate(
      post_allele_freqs_lower = 1,    # Add new columns with value 1
      post_allele_freqs_med = 1,
      post_allele_freqs_upper = 1,
      freq = 1
    ) %>%
    select(
      post_allele_freqs_lower,        # Reorder to place the new columns at the beginning
      post_allele_freqs_med,
      post_allele_freqs_upper,
      freq,
      everything()
    )

  allele_freq_summary <- rbind(allele_freq_summary, one_allele_loci)

  # Add removed loci to he_summary
  target_id_count <- moire_object$moire_data %>%
    select(!sample_id) %>%
    group_by(locus) %>%         # Group by 'locus' and 'allele'
    summarise(sample_total = n(), .groups = 'drop') %>% # Count the occurrences and drop grouping
    rename(target_id = locus)

  he_summary <- target_id_count %>%
    full_join(he_summary, by = "target_id")

  # Write summaries to files
  readr::write_tsv(coi_summary, coi_summary_o)
  readr::write_tsv(he_summary, he_summary_o, na = "0")
  readr::write_tsv(allele_freq_summary, allele_freq_summary_o)
  readr::write_tsv(relatedness_summary, relatedness_summary_o)
  readr::write_tsv(effective_coi_summary, effective_coi_summary_o)
}

# Main-----------------------------------------------------------------

# Create Moire object -------------------------------------------------
moire_object <- create_moire_input(arg$allele_table,
  arg$allow_relatedness,
  arg$burnin,
  arg$samples_per_chain,
  arg$verbose,
  arg$eps_pos_alpha,
  arg$eps_pos_beta,
  arg$eps_neg_alpha,
  arg$eps_neg_beta,
  arg$r_alpha,
  arg$r_beta,
  arg$mean_coi_shape,
  arg$mean_coi_scale,
  arg$max_eps_pos,
  arg$max_eps_neg,
  arg$record_latent_genotypes,
  arg$pt_chains,
  arg$pt_grad_lower,
  arg$pt_num_threads,
  arg$adapt_temp,
  arg$max_runtime
)

# Run Moire -------------------------------------------------------------------

moire_results <- run_moire(moire_object)

# Generate summaries
summarize_and_write_results(moire_object, moire_results, arg$coi_summary, arg$he_summary, arg$allele_freq_summary, arg$relatedness_summary, arg$effective_coi_summary)

