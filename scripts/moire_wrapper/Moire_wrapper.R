library(tibble)
library(dplyr)
library(moire)
library(tidyr)
library(optparse)
library(stringr)
library(validate)
library(checkmate)

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
      "Shape parameter for the mean coefficient of inbreeding prior.",
      "Default set to 0.1."
    )
  ),
  make_option(
    "--mean_coi_scale",
    type = "integer",
    default = 10,
    help = str_c(
      "Scale parameter for the mean coefficient of inbreeding prior.",
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
    "--num_chains",
    type = "integer",
    default = 1,
    help = str_c(
      "Number of chains to run in the MCMC process.",
      "Default set to 1."
    )
  ),
  make_option(
    "--num_cores",
    type = "integer",
    default = 1,
    help = str_c(
      "Number of cores to use for parallel processing.",
      "Default set to 1."
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
    "--pt_grad",
    type = "integer",
    default = 1,
    help = str_c(
      "Number of gradient evaluations to use in the parallel tempering process.",
      "Default set to 1."
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
    default = FALSE,
    help = str_c(
      "Logical flag to enable adaptive temperature in the MCMC process.",
      "Default set to FALSE."
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
    "--seed",
    type = "integer",
    default = 1,
    help = str_c(
      "Random number seed",
      "Default set to 1")
  )
)

arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table.tsv"
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
#' @param num_chains Numeric. Number of MCMC chains.
#' @param num_cores Numeric. Number of CPU cores to use for parallel processing.
#' @param pt_chains Numeric. Number of chains for parallel tempering.
#' @param pt_grad Numeric. Gradient step size for parallel tempering.
#' @param pt_num_threads Numeric. Number of threads for parallel tempering computations.
#' @param adapt_temp Logical. Whether to adapt the temperature during parallel tempering.
#' @param max_runtime Numeric. Maximum runtime allowed for the MCMC algorithm (in seconds).
#' @param seed Numeric or NULL. Random seed for reproducibility. If NULL, no seed is set.
#'
#' @return A list containing the Moire data and parameters, ready for downstream analysis.
#' @examples
#' \dontrun{
#' moire_input <- create_Moire_input(
#'   input_path = "data.csv", allow_relatedness = TRUE, burnin = 1000, 
#'   samples_per_chain = 5000, verbose = TRUE, eps_pos_alpha = 2, 
#'   eps_pos_beta = 5, eps_neg_alpha = 2, eps_neg_beta = 5, r_alpha = 1, 
#'   r_beta = 1, mean_coi_shape = 2, mean_coi_scale = 0.5, max_eps_pos = 0.1, 
#'   max_eps_neg = 0.1, record_latent_genotypes = FALSE, num_chains = 4, 
#'   num_cores = 2, pt_chains = 2, pt_grad = 0.01, pt_num_threads = 2, 
#'   adapt_temp = TRUE, max_runtime = 3600, seed = 42
#' )
#' }
#' @export

create_Moire_input <- function(input_path, allow_relatedness, burnin, 
                               samples_per_chain, verbose, eps_pos_alpha, 
                               eps_pos_beta, eps_neg_alpha, eps_neg_beta, 
                               r_alpha, r_beta, mean_coi_shape, mean_coi_scale, 
                               max_eps_pos, max_eps_neg, record_latent_genotypes, 
                               num_chains, num_cores, pt_chains, pt_grad, 
                               pt_num_threads, adapt_temp, max_runtime, seed) {
  print("Reading input data")
  input_data <- read.csv(input_path, na.strings = "NA", sep = "\t")
  moire_data <- input_data |>
    dplyr::select(specimen_id, target_id, seq) |> 
    dplyr::rename(sample_id = specimen_id, locus = target_id, allele = seq)

  print("Creating Moire object")
  # Create a list containing selected data and parameters
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
      num_chains = num_chains,
      num_cores = num_cores,
      pt_chains = pt_chains,
      pt_grad = pt_grad,
      pt_num_threads = pt_num_threads,
      adapt_temp = adapt_temp,
      max_runtime = max_runtime,
      seed = seed
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
  assert_numeric(moire_object$moire_parameters$num_chains, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$num_cores, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$pt_chains, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$pt_grad, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$pt_num_threads, any.missing = FALSE, len = 1)
  assert_logical(moire_object$moire_parameters$adapt_temp, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$max_runtime, any.missing = FALSE, len = 1)
  assert_numeric(moire_object$moire_parameters$seed, any.missing = FALSE, null.ok = TRUE, len = 1)  # Allow NULL seed
  
  
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

run_Moire <- function(moire_object) {

  moire_data = moire_object$moire_data |>
    moire::load_long_form_data()

  moire_parameters = moire_object$moire_parameters

  runtime <- system.time({
    set.seed(moire_parameters$seed)
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
      num_chains = moire_parameters$num_chains,
      num_cores = moire_parameters$num_cores,
      pt_chains = moire_parameters$pt_chains,
      pt_grad = moire_parameters$pt_grad,
      pt_num_threads = moire_parameters$pt_num_threads,
      adapt_temp = moire_parameters$adapt_temp,
      max_runtime = moire_parameters$max_runtime
    )
  })
  
  return(moire_res)
}

#' Summarize and Write MCMC Results
#'
#' This function processes MCMC analysis results to generate summaries of key statistics and writes the results to TSV files.
#'
#' @param mcmc_results A list containing the results of the MCMC analysis, structured for input to the `moire` package's summarization functions.
#'
#' @details
#' The function computes the following summaries:
#' \describe{
#'   \item{\code{coi_summary}}{Summarizes the complexity of infection (COI) for each specimen, with columns:
#'     \itemize{
#'       \item \code{specimen_id}: Identifier for each specimen.
#'       \item \code{coi}: Mean posterior COI.}}
#'   \item{\code{he_summary}}{Summarizes heterozygosity for each target, with columns:
#'     \itemize{
#'       \item \code{target_id}: Identifier for each genetic locus.
#'       \item \code{he}: Mean posterior heterozygosity. Will be NA if locus was filtered out by Moire (either because it is uninformative or because all but one of the alleles are very low frequency.}}
#'   \item{\code{allele_freq_summary}}{Summarizes allele frequencies for each target, with columns:
#'     \itemize{
#'       \item \code{target_id}: Identifier for each genetic locus.
#'       \item \code{freq}: Mean posterior allele frequency.}}
#'   \item{\code{relatedness_summary}}{Summarizes within-host relatedness for each specimen, with columns:
#'     \itemize{
#'       \item \code{specimen_id}: Identifier for each specimen.
#'       \item \code{within_host_rel}: Mean posterior within-host relatedness.}}
#'   \item{\code{effective_coi_summary}}{Summarizes the effective complexity of infection for each specimen, with columns:
#'     \itemize{
#'       \item \code{specimen_id}: Identifier for each specimen.
#'       \item \code{ecoi}: Mean posterior effective COI.}}
#' }
#'
#' The function writes the summaries to the following files:
#' \itemize{
#'   \item \code{"coi_summary.tsv"}
#'   \item \code{"he_summary.tsv"}
#'   \item \code{"allele_freq_summary.tsv"}
#'   \item \code{"relatedness_summary.tsv"}
#'   \item \code{"effective_coi_summary.tsv"}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' summarize_and_write_results(mcmc_results)
#' }
#'
#' @importFrom dplyr rename
#' @importFrom readr write_tsv
#' @importFrom moire summarize_coi summarize_he summarize_allele_freqs summarize_relatedness summarize_effective_coi
#' @export
summarize_and_write_results <- function(moire_object, mcmc_results) {
  
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
  missing_target_ids = unique(moire_object$moire_data$locus[
    !moire_object$moire_data$locus %in% allele_freq_summary$target_id])
  present_target_ids = unique(moire_object$moire_data$locus[
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
  
  allele_freq_summary = rbind(allele_freq_summary,
        one_allele_loci)
  
  # Add removed loci to he_summary
  target_id_count <- moire_object$moire_data %>% 
    select(!sample_id) %>%
    group_by(locus) %>%         # Group by 'locus' and 'allele'
    summarise(sample_total = n(), .groups = 'drop') %>% # Count the occurrences and drop grouping
    rename(target_id = locus)
  
  he_summary <- target_id_count %>%
    full_join(he_summary, by = "target_id")
  
  # Write summaries to files
  readr::write_tsv(coi_summary, "coi_summary.tsv")
  readr::write_tsv(he_summary, "he_summary.tsv")
  readr::write_tsv(allele_freq_summary, "allele_freq_summary.tsv")
  readr::write_tsv(relatedness_summary, "relatedness_summary.tsv")
  readr::write_tsv(effective_coi_summary, "effective_coi_summary.tsv")
}

# Main-----------------------------------------------------------------

# Create Moire object -------------------------------------------------
moire_object <- create_Moire_input(arg$allele_table,
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
  arg$num_chains,
  arg$num_cores,
  arg$pt_chains,
  arg$pt_grad,
  arg$pt_num_threads,
  arg$adapt_temp,
  arg$max_runtime,
  arg$seed
)

# Run Moire -------------------------------------------------------------------

moire_results <- run_Moire(moire_object)

# Generate summaries
summarize_and_write_results(moire_object, moire_results)

