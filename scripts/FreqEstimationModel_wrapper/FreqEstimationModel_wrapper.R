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
library(readr)
library(purrr)
library(variantstring)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--aa_calls", 
    help = str_c(
      "TSV containing amino acid calls, with the columns: specimen_id, ", 
      "target_id, gene_id, aa_position, ref_codon, ref_aa, codon, aa. Required."
    )
  ), 
  make_option(
    "--coi", 
    help = str_c(
      "TSV containing COI for each specimen, with the columns: specimen_id, ", 
      "coi, or numeric input of COI. Required."
    )
  ), 
  make_option(
    "--groups", 
    help = str_c(
      "TSV containing which regions should be analyzed for each group, with the columns:
      group_id, gene_id, aa_position. Required."
    )
  ),
  make_option("--seed", type = "integer", help = "Random number seed. Optional.", default=1), 
  make_option(
    "--mlaf_output", 
    help = str_c(
      "Output TSV to contain multilocus allele frequencies, with the columns: ", 
      "variant, freq, median_freq, CI_2.5, CI_97.5, prev, sample_total. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))


if(interactive()){
  arg$groups <- "example_loci_groups.tsv"
  arg$coi <- "example_coi_table.tsv"
  arg$seed <-  "1"
  arg$aa_calls <- "example_amino_acid_calls.tsv"
  arg$mlaf_output <- "output/FEM_output.tsv"
}

#' Returns average COI from COI tsv file path
#'
#' Takes in the COI table for the specimens and calculates the average COI,
#' or returns the COI value if given directly.
#'
#' @param coi_path Path of TSV file with COI for each specimen. It should 
#' have two columns, "specimen_id" and "coi." Alternatively, you can provide
#' the numeric value of the average COI. Parsing of which option is done 
#' automatically.
#' 
#' @return average COI across all samples in the file
calculate_avg_COI <- function(coi_path){
  coi <- as.numeric(coi_path)
  if(!is.na(coi)){
    return(as.numeric(coi_path))
  }
  else{ #if NA with casting the string to numeric, we were given a path
    COI_table <- read_tsv(coi_path, col_types=list(
      "specimen_id"=col_character(),
      "coi"=col_number()))
    rules <- validate::validator(
      ! is.na(specimen_id), 
      ! is.na(coi)
    )
    fails <- validate::confront(COI_table, rules, raise = "all") %>%
      validate::summary() %>%
      dplyr::filter(fails > 0)
    if (nrow(fails) > 0) {
      stop(
        "Input input_data failed one or more validation checks: ", 
        str_c(fails$expression, collapse = "\n"), 
        call. = FALSE
      )
    }
  }
  
  vals <- COI_table$coi
  avg_val <- mean(vals)
  return(avg_val)
}

#' Reads in group file
#'
#' Takes in the group file path and returns a dataframe object
#'
#' @param groups_path Path of TSV file with COI for each specimen. It should have
#' three columns: group_id, gene_id, aa_position
#' 
#' @return a group_table dataframe
read_groups <- function(groups_path){
  group_table <- read_tsv(groups_path, col_types=list(
    "group_id"=col_character(),
    "gene_id"=col_character(),
    "aa_position"=col_integer()
  ))
  return(group_table)
}

#' Returns dataframe with only biallelic SNPS
#'
#' Takes a dataframe in the long-data format and returns only the targets that
#' are biallelic in the population
#'
#' @param input_data dataframe object containing the columns specimen_id, target_id,
#' read_count, gene_id, aa_position, ref_codon, ref_aa, codon, aa
#' 
#' @return that same dataframe object but with only bi or mono-allelic targets
check_biallelic <- function(input_data){
  mutants <- input_data %>% filter(ref_aa != aa,) %>%
    distinct(unique_targets, aa, keep.all=T)
  bad_targets <- mutants %>% filter(unique_targets %in% duplicated(unique_targets)) #remove targets that have >1 non-reference call
  if(nrow(bad_targets)>0){
    warning("Not biallelic, dropped", bad_targets)
  }
  only_biallelic <- input_data %>% filter(!(unique_targets %in% bad_targets))
  alt_calls <- only_biallelic %>% 
    filter(ref_aa != aa) #dataframe containing alt + ref calls for biallelic SNPs
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
#' @param groups output of read_groups
#' @param group_id the string of the current group being analyzed
#' 
#' @return Matrix of frequencies for each multi-locus genotype.
create_FEM_input <- function(input_path, groups, group_id) {
  input_data <- read_tsv(input_path, col_types=list("specimen_id"=col_character(),
                                                    "aa_position" = col_integer(),
                                                    "target_id"=col_character(),
                                                    "gene_id"=col_character(),
                                                    "ref_codon"=col_character(),
                                                    "ref_aa"=col_character(),
                                                    "codon"=col_character(),
                                                    "aa"=col_character()))
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
  input_data$unique_targets <- paste(input_data$gene_id, input_data$aa_position, sep=":")
  groups <- groups[groups$group_id == group_id,]
  group_targets <- paste(groups$gene_id, groups$aa_position, sep=":")
  group_input_data <- input_data[input_data$unique_targets %in% group_targets,]
  input_data <- group_input_data
  data_list <- check_biallelic(input_data)
  input_data <- data_list[[1]]
  alt_alleles <- data_list[[2]]

  unique_targets <- unique(input_data$unique_targets)
  unique_sample_ids <- unique(input_data$specimen_id)
  
  sample_matrix <- matrix(99,  nrow=length(unique_sample_ids), ncol=length(unique_targets))
  #99 is the "no data" identifier
  colnames(sample_matrix) <- unique_targets
  rownames(sample_matrix) <- unique_sample_ids
 
  # populate the input matrix with 0 (only reference), 0.5 (het cal), or 1 (only non-reference)
  for(sample in input_data$specimen_id){
    cut_df <- input_data[input_data$specimen_id==sample,]
    for(unique_target in unique_targets){
      cut_cut_df <- cut_df[cut_df$unique_targets==unique_target,]
      nvals <- 99
      nvals <- length(unique(cut_cut_df$aa)) 
      if (nvals > 2){
        stop(
          "Too many alleles for FEM", 
          call. = FALSE
        )
      }
      if (nvals == 2){
        sample_matrix[sample, unique_target] <- 0.5 # "heterozygous" call
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
  num_group <- nrow(sample_matrix)
  return(list(
    sample_matrix,
    alt_alleles,
    num_group))
}
 
#' Run FEM
#'
#' Given the input_data path, groups file, group_id, and average COI, runs
#' FreqEstimationModel
#' 
#' @param input_data_path Path of TSV file with amino acid calls. It should 
#'   have character columns for specimen_id, target_id, gene_id, 
#'   ref_codon, ref_aa, codon, and aa. It should have an integer 
#'   aa_position column. There should be no explicit missing data.
#' @param groups output of read_groups
#' @param coi output of calculate_avg_COI
#' @param group single group being analyzed
#' 
#' @return list of 4 elements: plsf_table, runtime information, target_mapping, and
#' alternative alleles
run_FreqEstimationModel <- function(input_data_path, groups, COI, group) {
    sample_matrix_list <- create_FEM_input(input_data_path, groups, group)
    sample_matrix <- sample_matrix_list[[1]]
    alt_alleles <- sample_matrix_list[[2]]
    num_group <- sample_matrix_list[[3]]
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
        moi_hyperparameter <- COI 
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
        
        set.seed(as.numeric(arg$seed))

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
        pop_prev <- 1-(1-pop_freq)^median(mcmc_mois)
        inf_prev <- t(apply(mcmc_As, 2, function(x) colMeans(x > 0)))
        prev <- colMeans(inf_prev)
        pop_freq <- cbind(pop_freq, prev)
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
            alt_allele = alt_alleles,
            num_group = num_group
        )
    )
}
#' Turns a binary string into the STAVE format
#'
#' Takes a string outputted as "sequence" from FEM and given the mapping of targets
#' to positions in that string, returns a STAVE-readable format
#'
#' @param chars binary string representing each sequence
#' @param names output from run_FEM mapping positions in the string to targets
#' @param alt_alleles also an output from run_FEM containing ref + alt allele identities
#' 
#' @return string of the STAVE format for a given string
bin2STAVE <- function(chars, names, alt_alleles){
  char_ix <- 1
  chars_split <- strsplit(chars, "")[[1]]
  gene_list <- c()
  long_form <- tibble(
    gene = character(),
    pos = numeric(),
    n_aa = numeric(),
    het = logical(),
    phased = logical(),
    aa = character(),
    read_count = numeric()
  )
  
  for(char in chars_split){
    if(char==1){call <- "aa"}
    else {call <- "ref_aa"}
    alt_current <- alt_alleles[alt_alleles$unique_targets == names[char_ix],]
    gene_current <- str_split(names[char_ix], ":")[[1]][1] ###split appropriately
    pos_current <- str_split(names[char_ix], ":")[[1]][2] ###split appropriately
    gene_list <- append(gene_list, gene_current)
    alt_current <- alt_current[1, call]
    alt_current <- as.character(alt_current)
    char_ix <- char_ix + 1
    long_form <- long_form %>% 
      add_row(gene = gene_current,
              pos = as.numeric(pos_current),
              n_aa = NA,
              het = NA,
              phased = TRUE,
              aa = alt_current,
              read_count = NA,)
  }
  #single_locus_STAVE, into multi-locus stave
  #adding if calls are het or not in the given locus
  long_form <- long_form %>% 
    group_by(gene, pos) %>%
    mutate(n_aa = n())
    
  long_form <- long_form %>%
    mutate(het = if_else(n_aa>1, TRUE, FALSE))
  string_output <- long_to_variant(list(long_form))
  return(string_output)
}

#' Formats a single output from run_FEM
#'
#' Takes the output of run_FreqEstimationModel and formats to include the STAVE format
#'
#' @param pop_freq_list output list from run_FreqEstimationModel
#' 
#' @return formatted dataframe
format_single_group_output <- function(pop_freq_list){
  input_list <- pop_freq_list[[1]]
  names <- pop_freq_list[[3]]
  alt_alleles <- pop_freq_list[[4]]
  num_group <- pop_freq_list[[5]]
  
  input_list <- input_list %>% 
    mutate(variant=map_chr(sequence, bin2STAVE, names, alt_alleles),
                        sample_total=num_group) %>% 
    as_tibble()
  return(input_list)
}


groups <- read_groups(arg$groups)
COI <- calculate_avg_COI(arg$coi)
overall_output <- data.frame("sequence"=c(),	"freq"=c(),	"median_freq"=c(),
                             "CI_2.5"=c(),	"CI_97.5"=c())
#run FEM separately for each group
for(group in unique(groups$group_id)){
  #will be one output file
  fem_results <- run_FreqEstimationModel(arg$aa_calls, groups, COI, group)
  fem_plsf <- format_single_group_output(fem_results)
  fem_plsf$group_id <- group
  overall_output <- rbind(overall_output, fem_plsf)
}
#format and write output to disk
overall_output <- apply(overall_output,2,as.character)
overall_output_df <- data.frame(overall_output)
write_tsv(overall_output_df, arg$mlaf_output)



