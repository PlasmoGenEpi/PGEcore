library(tibble)
library(dplyr)
library(FreqEstimationModel)
library(foreach)
library(doMC)
library(iterators)
library(parallel)
library(rngtools)
library(validate)
library(magrittr)
library(stringr)
library(optparse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--aa_calls", 
    help = str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ", 
      "target_id, gene_id, aa_position, ref_codon, ref_aa, codon, aa"
    )
  ), 
  make_option(
    "--coi", 
    help = str_c(
      "TSV containing COI for each specimen, with the columns: specimen_id, ", 
      "coi"
    )
  ), 
  make_option("--seed", type = "integer", help = "Random number seed"), 
  make_option(
    "--mlaf", 
    help = str_c(
      "TSV containing multilocus allele frequencies, with the columns: ", 
      "variant, freq, total"
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))

#' Reformat amino acid calls into the form required by FEM
#'
#' Takes the path of TSV file of amino acid calls, reads it into a 
#' tibble, and reformats that tibble into the form needed by FEM.
#'
#' @param input_path Path of TSV file with amino acid calls. It should 
#'   have character columns for specimen_id, target_id, gene_id, 
#'   ref_codon, ref_aa, codon, and aa. It should have an integer 
#'   aa_position column. There should be no explicit missing data.
#' 
#' @return Matrix of frequencies for each multi-locus genotype.
create_FEM_input <- function(input_path) {
  input_data <- read.csv(input_path, na.strings = "NA")
  
  # Validate input format
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(target_id), 
    is.character(gene_id), 
    is.integer(aa_position), 
    is.character(ref_codon), 
    is.character(ref_aa), 
    is.character(codon), 
    is.character(aa), 
    ! is.na(specimen_id), 
    ! is.na(target_id), 
    ! is.na(gene_id), 
    ! is.na(aa_position), 
    ! is.na(ref_codon), 
    ! is.na(ref_aa), 
    ! is.na(codon), 
    ! is.na(aa)
  )
  fails <- validate::confront(input_data, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }
  
  input_data$unique_targets <- paste(input_data$target, input_data$position)
  unique_targets <- unique(input_data$unique_targets)
  unique_sample_ids <- unique(input_data$sample_id)
  
  sample_matrix <- matrix(99,  nrow=length(unique_sample_ids), ncol=length(unique_targets))
  colnames(sample_matrix) <- unique_targets
  rownames(sample_matrix) <- unique_sample_ids
  
  for(sample in input_data$sample_id){
    cut_df <- input_data[input_data$sample_id==sample,]
    for(unique_target in unique_targets){
      print(paste("Testing sample", sample, "at target", unique_target))
      cut_cut_df <- cut_df[cut_df$unique_targets==unique_target,]
      print(cut_cut_df)
      nvals <- 99
      nvals <- length(unique(cut_cut_df$codon)) ## DO WE CARE ABOUT SILENT MUTANTS?
      print(unique(cut_cut_df$AA))
      print(paste("nvals", nvals))
      if (nvals >= 2){
        sample_matrix[sample, unique_target] <- 0.5
      }
      if (nvals == 1){
        sample_matrix[sample, unique_target] <- 1
      }
    
    }
  }
  return(sample_matrix)
}
 

run_FreqEstimationModel <- function(input_data_path, output_dir) {
    # TODO: If missing data fill with 99
    tmp_file_path <- file.path(output_dir, "tmp_formatted_output_fem.txt")
    sample_matrix <- create_FEM_input(input_data_path, tmp_file_path)
    data_summary <- list()
    #data_summary$Data <- read.delim(tmp_file_path, row.names = 1) # Specify row.names = 1 to use the first column as row names
    #data_summary$Data <- as.matrix(data_summary$Data)
    data_summary$Data <- sample_matrix
    runtime <- system.time({
        thinning_interval <- 1 # Number of iterations per chain that are not saved
        no_traces_preburnin <- 10000 # For more traces but manageable pdf plots, don't exceed 10k and increase thinning interval instead
        no_mcmc_chains <- 3 # Number of MCMC chains to run
        parallel <- FALSE # Set to true if running code in parallel (this option isn't active yet)
        NGS <- FALSE # Set to true if data are in NGS format (this option isn't active yet)
        log_like_zero <- FALSE # QC check: sets log(p(yi|ai)) to zero (only impacts NGS)
        mcmc_variable_list <- list(
            no_mcmc_chains = no_mcmc_chains,
            no_traces_preburnin = no_traces_preburnin,
            thinning_interval = thinning_interval,
            NGS = NGS,
            log_like_zero = log_like_zero
        )
        if (!NGS | log_like_zero) { # Set to true to keep partial observations rather than discard them
            augment_missing_data <- TRUE
        } else {
            augment_missing_data <- FALSE
        }
        if (log_like_zero) { # When log_like_zero == TRUE, code should return the prior
            moi_prior <- "Uniform" # Use a prior that is easy to eye-ball
        } else {
            moi_prior <- "Poisson" # Choose between 'Uniform', 'Poisson', 'Geometric' or 'nBinomial'
        }
        moi_max <- 8 # Maximum MOI regarded as possible by the model (I haven't tested beyond 20)
        moi_hyperparameter <- 3 # Specify population-average MOI (parameter of the prior on the MOI)
        moi_size_hyperparameter <- 0.5 # Only applies if moi_prior == 'nBinomial' (hyperparameter for the prior on the MOI if the prior is the negative Binomial)
        moi_prior_min2 <- NULL # Specify the lower bound for the MOI per individual
        moi_initial <- NULL # If null, the initial vector of MOIs is set internally, otherwise set to input moi_initial
        moi_list <- list(
            moi_hyperparameter = moi_hyperparameter,
            moi_size_hyperparameter = moi_size_hyperparameter,
            moi_prior = moi_prior,
            moi_max = moi_max,
            moi_prior_min2 = moi_prior_min2,
            moi_initial = moi_initial
        )
        processed_data_list <- preprocess_data(
            data_summary,
            log_like_zero,
            NGS,
            augment_missing_data,
            moi_prior_min2
        )
        frequency_hyperparameter <- rep(1, processed_data_list$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency vector
        frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
        frequency_list <- list(
            frequency_hyperparameter = frequency_hyperparameter,
            frequency_initial = frequency_initial
        )
        seed <- 1 # For reproducibility
        set.seed(seed)

        # Run MCMC - works up to here
        results <- mcmc_sampling_parallel(processed_data_list,
            moi_list,
            frequency_list,
            mcmc_variable_list,
            cores_max = 3
        )

        # Generate numerical approximations of posteriors by removing burnin
        burnin <- 1:(0.5 * mcmc_variable_list$no_traces_preburnin) # remove first half
        if (mcmc_variable_list$no_mcmc_chains > 1) {
            alply_genotype_freq_store_chains_burnin <- plyr::alply(results$genotype_freq_store_chains[-burnin, , ], 3)
        } else {
            alply_genotype_freq_store_chains_burnin <- results$genotype_freq_store_chains[-burnin, , ]
        }

        # Use mcmc.list to represent frequency chains across parallel runs (see ?coda::mcmc.list)
        mcmc_frequency_chains <- coda::mcmc.list(lapply(alply_genotype_freq_store_chains_burnin, # 3 for splitting by third dimension (nothing to do with no. of chains)
            coda::mcmc,
            start = (max(burnin) + 1) * mcmc_variable_list$thinning_interval,
            end = mcmc_variable_list$no_traces_preburnin * mcmc_variable_list$thinning_interval,
            thin = mcmc_variable_list$thinning_interval
        ))

        mcmc_As <- abind::abind(plyr::alply(results$genotype_count_store_chains[-burnin, , , ], 4), along = 1) # Haplotype counts for all chains excluding burnin
        mcmc_mois <- apply(mcmc_As, c(1, 2), sum)

        # Extract population-level haplotype frequency estimates
        pop_freq <- cbind(
            fem_mean_frequency = summary(mcmc_frequency_chains)$statistics[, "Mean"],
            fem_median_frequency = summary(mcmc_frequency_chains)$quantiles[, 3],
            "fem_CI2.5%" = summary(mcmc_frequency_chains)$quantiles[, 1],
            "fem_CI97.5%" = CI_upper <- summary(mcmc_frequency_chains)$quantiles[, 5]
        )
    })
    # Convert row names to a new column named "sequence"
    sequence_column <- data.frame(sequence = rownames(pop_freq), stringsAsFactors = FALSE)
    pop_freq <- cbind(sequence_column, pop_freq)
    rownames(pop_freq) <- NULL
    file.remove(tmp_file_path)
    return(
        list(
            plsf_table = pop_freq,
            runtime = runtime
        )
    )
}

output_dir <- "output"
input_dir <- "test_dataset.csv"
fem_results <- run_FreqEstimationModel(input_dir, output_dir)
fem_plsf <- fem_results$plsf_table
write.csv(fem_plsf, file.path(output_dir, "plsf_table_STANDARDIZED.csv"))
