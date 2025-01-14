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
  make_option(
    "--groups", 
    help = str_c(
      "TSV containing which regions should be analyzed for each group, with the columns:
      group_id, gene_id, aa_position"
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
arg <- list(groups = "example_loci_groups.tsv",
            coi = "example_coi_table.tsv",
            aa_calls = "example_amino_acid_calls.tsv",
            seed = 1)


# be clear takes input path + returns number
calculate_COI <- function(coi_path){
  COI_table <- read.csv(coi_path, sep="\t")
  vals <- COI_table$coi
  avg_val <- mean(vals)
  return(avg_val)
}

read_groups <- function(groups_path){
  group_table <- read.csv(groups_path, sep="\t")
  return(group_table)
}


check_biallelic <- function(input_data){
  mutants <- input_data[input_data$ref_aa != input_data$aa,]
  mutants_ignoring_sample <- distinct(mutants, unique_targets, aa, keep.all=T)
  mutants <- mutants_ignoring_sample
  bad_targets <- mutants$unique_targets[duplicated(mutants$unique_targets)]
  print(paste('Dropped', bad_targets))
  only_biallelic <- input_data[!(input_data$unique_targets %in% bad_targets),]
  alt_calls <- only_biallelic[only_biallelic$ref_aa != only_biallelic$aa,]
  return(list(
    only_biallelic,
    alt_calls))
}

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
create_FEM_input <- function(input_path, group_id) {
  input_data <- read.csv(input_path, na.strings = "NA", sep="\t")
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
  
  input_data$unique_targets <- paste(input_data$gene_id, input_data$aa_position, sep=";")
  groups <- groups[groups$group_id == group_id,]
  group_targets <- paste(groups$gene_id, groups$aa_position, sep=";")
  group_input_data <- input_data[input_data$unique_targets %in% group_targets,]
  input_data <- group_input_data
  data_list <- check_biallelic(input_data)
  input_data <- data_list[[1]]
  alt_alleles <- data_list[[2]]

  unique_targets <- unique(input_data$unique_targets)
  unique_sample_ids <- unique(input_data$specimen_id)
  
  sample_matrix <- matrix(99,  nrow=length(unique_sample_ids), ncol=length(unique_targets))
  colnames(sample_matrix) <- unique_targets
  rownames(sample_matrix) <- unique_sample_ids
  
 
  
  for(sample in input_data$specimen_id){
    cut_df <- input_data[input_data$specimen_id==sample,]
    for(unique_target in unique_targets){
      cut_cut_df <- cut_df[cut_df$unique_targets==unique_target,]
      nvals <- 99
      nvals <- length(unique(cut_cut_df$codon)) ## DO WE CARE ABOUT SILENT MUTANTS?
      if (nvals >= 2){
        sample_matrix[sample, unique_target] <- 0.5
      }
      if (nvals == 1){
        if (cut_cut_df[1, "ref_aa"]==cut_cut_df[1, "aa"]){
          val <- 0
        }
        else {val <- 1}
        sample_matrix[sample, unique_target] <- val
      }
    
    }
  }
  return(list(
    sample_matrix,
    alt_alleles))
}
 

run_FreqEstimationModel <- function(input_data_path, group, COI, output_dir, seed) {
    sample_matrix_list <- create_FEM_input(input_data_path, group)
    sample_matrix <- sample_matrix_list[[1]]
    alt_alleles <- sample_matrix_list[[2]]
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
        #### ALFRED 
        moi_hyperparameter <- COI # Specify population-average MOI (parameter of the prior on the MOI)
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
        
        test_process <- preprocess_data(data_summary, FALSE, FALSE, TRUE, FALSE)
        
        frequency_hyperparameter <- rep(1, processed_data_list$no_haplotypes) # The Parameter vector for the Dirichlet prior on the frequency vector
        frequency_initial <- NULL # If null, external_frequency_initial set internally, otherwise set to input external_frequency_initial
        frequency_list <- list(
            frequency_hyperparameter = frequency_hyperparameter,
            frequency_initial = frequency_initial
        )
        seed <- seed
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
            freq = summary(mcmc_frequency_chains)$statistics[, "Mean"],
            median_freq = summary(mcmc_frequency_chains)$quantiles[, 3],
            "CI_2.5" = summary(mcmc_frequency_chains)$quantiles[, 1],
            "CI_97.5" = CI_upper <- summary(mcmc_frequency_chains)$quantiles[, 5]
        )
    })
    # Convert row names to a new column named "sequence"
    sequence_column <- data.frame(sequence = rownames(pop_freq), stringsAsFactors = FALSE)
    pop_freq <- cbind(sequence_column, pop_freq)
    rownames(pop_freq) <- NULL
    return(
        list(
            plsf_table = pop_freq,
            runtime = runtime,
            names = processed_data_list[["markerID"]],
            alt_allele = alt_alleles
        )
    )
}

bin2STAVE <- function(chars, names, alt_alleles){
  name <- ""
  char_ix <- 1
  chars_split <- strsplit(chars, "")[[1]]
  for(char in chars_split){
    if(char==1){call <- "aa"}
    else {call <- "ref_aa"}
      #gene;position;aa:gene;position;aa
    formatted_name <- gsub(" ", ";", name)
    alt_current <- alt_alleles[alt_alleles$unique_targets == names[char_ix],]
    alt_current <- alt_current[1, call]
    name <- paste(formatted_name, paste(names[char_ix], alt_current, sep=";"), sep=":")
  
    char_ix <- char_ix + 1
  }
  name <- sub('.', '', name)
  return(name)
}
format_output <- function(pop_freq_list){
  input_list <- pop_freq_list[[1]]
  names <- pop_freq_list[[3]]
  alt_alleles <- pop_freq_list[[4]]
  
  input_list$variant <- lapply(input_list$sequence, bin2STAVE, names, alt_alleles)
  return(input_list)
}

COI <- calculate_COI(arg$coi)
group <- read_groups(arg$group)
output_dir <- "output"
input_dir <- arg$aa_calls
seed <- arg$seed

#TODO put into function
#TODO add "total" column
overall_output <- data.frame("sequence"=c(),	"freq"=c(),	"median_freq"=c(),
                             "CI_2.5"=c(),	"CI_97.5"=c())
for(group in unique(groups$group_id)){
  #will be one output file
  fem_results <- run_FreqEstimationModel(input_dir, group, COI, output_dir, seed)
  fem_plsf <- format_output(fem_results)
  fem_plsf$group_id <- group
  overall_output <- rbind(overall_output, fem_plsf)
}
overall_output <- apply(overall_output,2,as.character)
overall_output <- overall_output[, 2:7]
write.csv(unlist(overall_output), file.path(output_dir, "plsf_table_grouped.csv"), row.names=FALSE)

#TODO add more docstrings
#TODO add in-line code as needed