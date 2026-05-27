#!/usr/bin/env Rscript

# Estimate IBD-based relatedness with Dcifer

# Load required libraries ----------------------------------------------
library(dcifer, warn.conflicts = F)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(parallel, warn.conflicts = F)
library(foreach, warn.conflicts = F)
library(iterators, warn.conflicts = F)
library(doParallel, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(optparse, warn.conflicts = F)
library(parallelly, warn.conflicts = F)
library(purrr, warn.conflicts = F)
library(readr, warn.conflicts = F)
library(stringr, warn.conflicts = F)
library(tibble, warn.conflicts = F)
library(tidyr, warn.conflicts = F)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing alleles, with columns identifying specimens, ", 
      "target names, and target values. The names of these columns are given ", 
      "by the --specimen_name_col, --target_name_col, and --target_value_col ", 
      "arguments, respectively. Required."
    )
  ), 
  make_option(
    "--coi_table", 
    help = 
      str_c(
        "TSV containing specimen COIs, with a coi column and a column with ", 
        "specimen IDs named according to --specimen_name_col. Optional."
      )
  ), 
  make_option(
    "--allele_freq_table", 
    help = 
      str_c(
        "TSV containing single locus allele frequencies, with a column with ", 
        "target names, named according to --target_name_col, a column with ", 
        "target values, named according to --target_value_col, and freq and ", 
        "total columns. Optional."
      )
  ), 
  make_option(
    "--specimen_name_col", 
    default = "specimen_name", 
    help = "String giving the name of the specimen ID column"
  ), 
  make_option(
    "--target_name_col", 
    default = "target_name", 
    help = 
      "String giving the name of the target name column (e.g., the locus name)"
  ), 
  make_option(
    "--target_value_col", 
    default = "seq", 
    help = 
      str_c(
        "String giving the name of the target value column (e.g., the allele ", 
        "call)"
      )
  ), 
  make_option(
    "--specimen_metadata", 
    help = 
      str_c(
        "TSV containing specimen metadata. It should have a specimen name ", 
        "column, named according to --specimen_name_col, and a population ", 
        "name column. named according to --pop_name_col. Optional."
      )
  ), 
  make_option(
    "--pop_name_col", 
    help = str_c(
      "Column in specimen metadata to use for calculating allele ", 
      "frequencies. If this is provided, allele frequencies will be ", 
      "calculated separately for each population and pair of populations ", 
      "(for between population comparisons). Optional."
    )
  ), 
  make_option(
    "--rnull", 
    type = "double", 
    default = 0, 
    help = "Relatedness value to use as null for hypothesis testing. Optional."
  ), 
  make_option(
    "--alpha", 
    type = "double", 
    default = 0.05, 
    help = "alpha value to use for hypothesis testing. Optional."
  ), 
  make_option(
    "--use_estm", 
    action = "store_true", 
    default = FALSE, 
    help = "Use `dcifer::ibdEstM` in place of `dcifer::ibdPair`"
  ), 
  make_option(
    "--threads", 
    type = "integer", 
    default = 1, 
    help = "Number of threads to use. Optional."
  ), 
  make_option(
    "--seed", 
    type = "integer", 
    default = 1, 
    help = "Random number seed. Optional."
  ), 
  make_option(
    c("-v", "--verbose"), 
    action = "store_true", 
    default = FALSE, 
    help = "Print detailed process output"
  ), 
  make_option(
    "--btwn_host_rel_output", 
    help = str_c(
      "Path of TSV file to contain relatedness results, with the columns: ", 
      "specimen_name, btwn_host_rel. Depending on parameters, it may also ", 
      "contain columns for p_value, CI_lower, CI_upper, and/or strain_pair. ", 
      "Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example2_allele_table.tsv"
  arg$specimen_metadata <- "../../data/example2_specimen_metadata.tsv"
  arg$pop_name_col <- "population"
  arg$btwn_host_rel_output <- "../../btwn_host_rel.tsv"
  arg$threads <- 2
  arg$seed <- 1
  arg$verbose <- FALSE
}

#' Read allele table into a tibble
#'
#' Read the allele table TSV into a tibble
#'
#' @param allele_table_path Path to the TSV file containing the allele 
#'   table.
#' @param specimen_name_col String giving the name of the specimen ID 
#'   column.
#' @param target_name_col String giving the name of the target name column.
#' @param target_value_col String giving the name of the column 
#'   containing target values (i.e., the genotypes).
#'
#' @return A tibble containing columns for specimen_name, target_name, and 
#'   target_value.
create_allele_table_input <- function(
                                      allele_table_path, 
                                      specimen_name_col = "specimen_name", 
                                      target_name_col = "target_name", 
                                      target_value_col = "seq") {

  # Read in table
  allele_table <- read_tsv(
      allele_table_path, 
      col_types = cols(.default = col_character()), 
      progress = FALSE
    ) %>%
    select(all_of(c(specimen_name_col, target_name_col, target_value_col))) %>%
    # Standardize names
    rename(
      specimen_name = all_of(specimen_name_col), 
      target_name = all_of(target_name_col), 
      target_value = all_of(target_value_col)
    )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_name), 
    is.character(target_name), 
    is.character(target_value), 
    ! is.na(specimen_name), 
    ! is.na(target_name), 
    ! is.na(target_value)
  )
  fails <- validate::confront(allele_table, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  allele_table

}

#' Read COI table into format needed by Dcifer
#'
#' Read the TSV of specimen COIs into a tibble, join with specimen IDs 
#' from the allele list to ensure order matches, and return as a vector.
#'
#' @param coi_path Path to TSV containing specimen COIs. It should 
#'   have a character column with specimen names, with name given by 
#'   \code{specimen_name_col}, and an integer coi column.
#' @param allele_list The list format output by `dcifer::readDat` and 
#'   `dcifer::formatDat`.
#' @inheritParams create_allele_table_input
#'
#' @return Vector of COI values, one for each sample.
create_coi_input <- function(
                             coi_path, 
                             allele_list, 
                             specimen_name_col = "specimen_name") {

  # Read input table
  coi <- read_tsv(
      coi_path, 
      col_types = cols(.default = col_character(), coi = col_integer()), 
      progress = FALSE
    ) %>%
    select(all_of(specimen_name_col), coi) %>%
    rename(specimen_name = all_of(specimen_name_col))

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_name), 
    is.integer(coi), 
    ! is.na(specimen_name), 
    ! is.na(coi)
  )
  fails <- validate::confront(coi, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  # Check that all specimen IDs in the allele table are in the COI 
  # input, and vice versa
  specimen_inallele_notincoi <- setdiff(names(allele_list), coi$specimen_name)
  if (length(specimen_inallele_notincoi) > 0) {
    stop(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_inallele_notincoi, collapse = " "), 
      call. = FALSE
    )
  }
  specimen_incoi_notinallele <- setdiff(coi$specimen_name, names(allele_list))
  if (length(specimen_incoi_notinallele) > 0) {
    warning(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_incoi_notinallele, collapse = " "), 
      call. = FALSE
    )
  }

  # Join COI to allele list to ensure order matches
  coi <- tibble(specimen_name = names(allele_list)) %>%
    left_join(coi, by = "specimen_name")

  # Create and return named vector of COI
  setNames(coi$coi, coi$specimen_name)

}

#' Read allele frequency into a list format
#'
#' Read a TSV of allele frequencies into a tibble, then reformat into 
#' the list format taken by 'dcifer::ibdDat` and `dcifer::ibdPair`.
#'
#' @param allele_freq_path Path to the TSV file containing the allele 
#'   frequencies. There should be character columns for target names 
#'   and values, with names given by \code{target_name_col} and 
#'   \code{target_value_col}, respectively, and a double freq column.
#' @inheritParams create_allele_table_input, create_coi_input
#'
#' @return A named list with target_name as the names. The values are a 
#'   named vector with target_value as the name and freq as the value. This is 
#'   format taken by `dcifer::ibdDat` and `dcifer::ibdPair`.
create_allele_freq_input <- function(allele_freq_path, 
                                     allele_list, 
                                     target_name_col = "target_name", 
                                     target_value_col = "seq") {

  # Read in table
  allele_freqs <- read_tsv(
      allele_freq_path, 
      col_types = cols(
        .default = col_character(), 
        freq = col_double()
      ), 
      progress = FALSE
    ) %>%
    select(all_of(c(target_name_col, target_value_col)), freq) %>%
    # Standardize names
    rename(
      target_name = all_of(target_name_col), 
      target_value = all_of(target_value_col)
    )

  # Validate fields
  rules <- validate::validator(
    is.character(target_name), 
    is.character(target_value), 
    is.double(freq), 
    ! is.na(target_name), 
    ! is.na(target_value), 
    ! is.na(freq)
  )
  fails <- validate::confront(allele_freqs, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  # Check that all alleles in the allele frequency table are in the 
  # allele list, and vice versa
  allele_freq_alleles <- allele_freqs %>%
    unite(allele, target_name, target_value, sep = ":") %$%
    allele
  allele_table_alleles <- list_c(allele_list) %>%
    tibble(target_name = names(.), target_values = .) %>%
    mutate(target_values = map(target_values, names)) %>%
    unnest(target_values) %>%
    unite(alleles, target_name, target_values, sep = ":") %$%
    alleles
  specimen_intable_notinfreq <- setdiff(
    allele_table_alleles, 
    allele_freq_alleles
  )
  if (length(specimen_intable_notinfreq) > 0) {
    stop(
      "The following alleles are in the allele table and not in the provided ", 
      "allele frequencies: ", 
      str_c(specimen_intable_notinfreq, collapse = " "), 
      call. = FALSE
    )
  }
  specimen_infreq_notintable <- setdiff(
    allele_freq_alleles, 
    allele_table_alleles
  )
  if (length(specimen_infreq_notintable) > 0) {
    stop(
      "The following alleles are in the provided allele frequencies and not ", 
      "in the allele table: ", 
      str_c(specimen_infreq_notintable, collapse = " "), 
      call. = FALSE
    )
  }

  # Convert to Dcifer format and return
  dcifer::formatAfreq(
    allele_freqs, 
    lvar = "target_name", 
    avar = "target_value", 
    fvar = "freq"
  )

}

#' Read specimen metadata into a tibble
#'
#' Read the specimen metadata TSV into a tibble
#'
#' @param specimen_metadata_path Path to the TSV file containing the 
#'   specimen metadata.
#' @param coi Vector of COI values with specimen_name names.
#' @param pop_name_col String giving the name of the population name 
#'   column.
#' @inheritParams create_allele_table_input
#'
#' @return A tibble containing columns for specimen_name and population.
create_specimen_metadata_input <- function(
                                      specimen_metadata_path, 
                                      coi, 
                                      specimen_name_col = "specimen_name", 
                                      pop_name_col = "population") {

  # Read in table
  specimen_metadata <- read_tsv(
      specimen_metadata_path, 
      col_types = cols(.default = col_character()), 
      progress = FALSE
    ) %>%
    select(all_of(c(specimen_name_col, pop_name_col))) %>%
    # Standardize names
    rename(
      specimen_name = all_of(specimen_name_col), 
      population = all_of(pop_name_col)
    )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_name), 
    is.character(population), 
    ! is.na(specimen_name), 
    ! is.na(population)
  )
  fails <- validate::confront(specimen_metadata, rules, raise = "all") %>%
    validate::summary() %>%
    dplyr::filter(fails > 0)
  if (nrow(fails) > 0) {
    stop(
      "Input input_data failed one or more validation checks: ", 
      str_c(fails$expression, collapse = "\n"), 
      call. = FALSE
    )
  }

  # Make sure all specimens have metadata
  specimens <- names(coi)
  specimens_notinmeta <- setdiff(specimens, specimen_metadata$specimen_name)
  if (length(specimens_notinmeta) > 0) {
    stop(
      "Some specimen IDs are missing from metadata: ", 
      str_c(specimens_notinmeta, collapse = ","), 
      call. = FALSE
    )
  }

  specimen_metadata

}

#' Run Dcifer with population-specific allele frequencies
#'
#' Read in the allele table and, using the provided population 
#' information and COIs, create sample allele calls, allele frequency, 
#' and COI inputs for Dcifer. Run Dcifer for each population or pair of 
#' populations, and bind the results into one tibble.
#'
#' @param specimen_metadata Tibble of specimen metadata. It should have 
#'   a character specimen_name column and a character population column.
#' @param use_estm Boolean indicating whether `dcifer::ibdEstM` should 
#'   be used in place of `dcifer::ibdPair`.
#' @inheritParams create_allele_table_input, create_specimen_metadata_input
#' @param ... Parameters to be passed to `run_ibdpair`.
#'
#' @return A tibble with a character sample_a column, a character 
#'   sample_b column, and a double estimate column. If `pval = TRUE`, 
#'   there will also be a double p_value column. If `confint = TRUE`, 
#'   there will also be double CI_lower and CI_upper columns.
#'
#' @details In the case of population pairs, allele frequencies are the 
#'   mean of the allele frequencies in each separate population.
run_dcifer_bypop <- function(
                             allele_table_path, 
                             coi, 
                             specimen_metadata, 
                             specimen_name_col = "specimen_name", 
                             target_name_col = "target_name", 
                             target_value_col = "seq", 
                             use_estm = FALSE, 
                             ...) {


  # Function to compute mean allele frequencies from two lists of 
  # allele frequencies
  allele_freq_mean <- function(af1, af2) {
    mean_af <- list()
    loci <- unique(c(names(af1), names(af2)))
    for (l in loci) {
      # Get allele names
      alleles <- unique(c(names(af1[[l]]), names(af2[[l]])))
      locus_mean_freqs <- numeric()
      for (a in alleles) {
        # Get allele frequencies for this allele
        af1_af2 <- c(af1[[l]][a], af2[[l]][a])
        # Replace missing allele frequencies with 0
        af1_af2[is.na(af1_af2)] <- 0
        locus_mean_freqs[a] <- mean(af1_af2)
      }
      mean_af[[l]] <- locus_mean_freqs
    }
    mean_af
  }
  
  # Read in allele table
  allele_table <- create_allele_table_input(
    allele_table_path, 
    specimen_name_col = arg$specimen_name_col, 
    target_name_col = arg$target_name_col, 
    target_value_col = arg$target_value_col
  )

  # Join to allele table
  alleles_w_specimen_meta <- allele_table %>%
    left_join(specimen_metadata, by = "specimen_name")

  # Prepare inputs and run Dcifer for each population separately
  pops <- unique(alleles_w_specimen_meta$population)
  # Save allele frequencies here for use with population combinations 
  # below
  allele_freq_lists <- list()
  dcifer_res <- list()
  for (pop_oi in pops) {
    alleles_filtered <- alleles_w_specimen_meta %>%
      filter(population == pop_oi)
    # The way R indexing works ensures coi will have the same specimen 
    # order as in alleles_filtered
    pop_coi <- coi[unique(alleles_filtered$specimen_name)]
    pop_alleles <- dcifer::formatDat(
        alleles_filtered, 
        svar = "specimen_name", 
        lvar = "target_name", 
        avar = "target_value"
      )
    allele_freq_lists[[pop_oi]] <- dcifer::calcAfreq(
      pop_alleles, 
      pop_coi, 
      tol = 1e-5
    )
    if (arg$use_estm) {
      dcifer_res[[pop_oi]] <- run_ibdestm(
          pop_alleles, 
          pop_coi, 
          allele_freq_lists[[pop_oi]], 
          ...
        )
    } else {
      dcifer_res[[pop_oi]] <- run_ibdpair(
          pop_alleles, 
          pop_coi, 
          allele_freq_lists[[pop_oi]], 
          ...
        )
    }
  }

  # Prepare inputs and run Dcifer for each pair of populations
  if (length(pops) > 1) {
    pop_combos <- combn(pops, 2)
    # Iterate through population combinations
    for (cb in 1:dim(pop_combos)[2]) {
      pop_combo <- pop_combos[,cb]
      # Sort population names before merging into one identifier, to 
      # enable consistent reconstruction
      pop_oi <- pop_combo %>%
        sort() %>%
        str_c(collapse = "_")
      # Compute mean of all allele frequencies between these two 
      # populations
      pop_combo_af <- allele_freq_mean(
        allele_freq_lists[[pop_combo[1]]], 
        allele_freq_lists[[pop_combo[2]]]
      )
      pop_combo_alleles <- alleles_w_specimen_meta %>%
        filter(population %in% pop_combo) %>%
        dcifer::formatDat(
          svar = "specimen_name", 
          lvar = "target_name", 
          avar = "target_value"
        ) %>%
        # Reorder allele list to match locus and allele order of the 
        # allele frequency list
        dcifer::matchAfreq(pop_combo_af) %$%
        dsmp
      # Use names from allele list to ensure order of specimens matches
      pop_combo_coi <- coi[names(pop_combo_alleles)]
      # Build list of sample pairs here and pass to Dcifer. We only want 
      # to run it for pairs between the two populations, not within.
      samples_pop_a <- alleles_w_specimen_meta %>%
        filter(population == pop_combo[1]) %$%
        unique(specimen_name)
      samples_pop_b <- alleles_w_specimen_meta %>%
        filter(population == pop_combo[2]) %$%
        unique(specimen_name)
      sample_pairs <- expand.grid(samples_pop_a, samples_pop_b) %>%
        as_tibble() %>%
        rename(sample_a = Var1, sample_b = Var2) %>%
        # expand.grid() returns a factor, and indexing with factors 
        # produces inconsistent results depending on whether a list or 
        # vector is being indexed. Therefore, it is important to convert 
        # to a character vector before proceeding.
        mutate(
          sample_a = as.character(sample_a), 
          sample_b = as.character(sample_b)
        )
      if (arg$use_estm) {
        dcifer_res[[pop_oi]] <- run_ibdestm(
            pop_combo_alleles, 
            pop_combo_coi, 
            pop_combo_af, 
            sample_pairs = sample_pairs, 
            ...
          )
      } else {
        dcifer_res[[pop_oi]] <- run_ibdpair(
            pop_combo_alleles, 
            pop_combo_coi, 
            pop_combo_af, 
            sample_pairs = sample_pairs, 
            ...
          )
      }
    }
  }

  # Bind and return results from each population and population pair
  reduce(dcifer_res, bind_rows) %>%
    return()

}

#' Run `ibdPair` in parallel
#'
#' This function uses `dcifer::ibdPair` to estimate relatedness among 
#' all sample pairs in `dsmp`. It essentially replicates the 
#' functionality of `dcifer::ibdDat` while allowing for parallel 
#' execution.
#'
#' @inheritParams run_dcifer_bypop
#' @inheritParams dcifer::ibdDat
#' @param sample_pairs Tibble with two character columns, sample_a and 
#'   sample_b, specifying the sample combinations for which relatedness 
#'   should be estimated. If not provided, all possible sample pairs 
#'   will be run.
#' @param total_cores Integer specifying the number of cores to use.
#' @param verbose If TRUE, output from each parallel process will be 
#'   printed.
#'
#' @return A tibble with a character sample_a column, a character 
#'   sample_b column, and a double estimate column. If `pval = TRUE`, 
#'   there will also be a double p_value column. If `confint = TRUE`, 
#'   there will also be double CI_lower and CI_upper columns.
#'
#' @details
#' If `rnull` is 0 or 1, Dcifer performs a one-sided hypothesis test, 
#' but it does a two-sided test if it is between 0 and 1. Therefore, if 
#' the scientific hypothesis of interest is, e.g., *r* > 0.25, it is 
#' necessary to divide the *p*-values returned by this function by two, 
#' and set *p*-values for *r* estimates below `rnull` to some 
#' arbitrarily high number if using something like the Benjamini-
#' Hochberg correction for multiple testing. These corrections are left 
#' to the user - they are NOT done in this function.
#'
#' This function was originally written by Max Murphy and 
#' lightly modified by Alfred Hubbard.
run_ibdpair <- function(
                            dsmp, 
                            coi, 
                            afreq, 
                            sample_pairs = NULL, 
                            pval = TRUE, 
                            confint = FALSE, 
                            rnull = 0, 
                            alpha = 0.05, 
                            nr = 1000, 
                            reval = NULL, 
                            total_cores = NULL, 
                            verbose = FALSE) {

    if (confint) {
        mnewton <- FALSE
        tol <- NULL
    } else {
        mnewton <- TRUE
        tol <- 1 / nr
    }
    if (!mnewton) {
        if (!inherits(reval, "matrix")) {
            reval <- generateReval(M = 1, rval = reval, nr = nr)
        }
        neval <- ncol(reval)
        logr <- logReval(reval, M = 1, neval = neval)
    } else {
        neval <- NULL
    }
    inull <- if (mnewton || !pval) {
        NULL
    } else {
        which.min(abs(reval - rnull))
    }
    afreq <- lapply(afreq, log)
    nloc <- length(afreq)

    # Determine side parameter based on rnull
    if (identical(rnull, 0L)) {
      side <- "right"
    } else if (identical(rnull, 1L)) {
      side <- "left"
    } else {
      side <- "two-sided"
    }

    if (is.null(sample_pairs)) {
      sample_pairs_matrix <- combn(names(dsmp), 2)
      sample_pairs <- tibble(
        sample_a = sample_pairs_matrix[1,], 
        sample_b = sample_pairs_matrix[2,]
      )
    }

    if (is.null(total_cores)) {
        total_cores <- parallelly::availableCores() - 1
    }
    
    if (is.null(getDefaultCluster())) {
        if (verbose) {
          cl <- makeCluster(total_cores, outfile = "")
        } else {
          cl <- makeCluster(total_cores)
        }
        setDefaultCluster(cl)
        registerDoParallel(cl)
    } else {
        cl <- getDefaultCluster()
    }

    res <- foreach(
            i = 1:total_cores, 
            .combine = rbind, 
            .packages = c("dcifer", "foreach", "iterators")
          ) %dopar% {
        total_pairs <- nrow(sample_pairs)
        # Use floor to avoid off-by-one error when these equations 
        # don't yield a whole number. All pairs will still be included 
        # because end_idx will be a whole number when i = total_cores.
        begin_idx <- floor(((i - 1) * total_pairs / total_cores) + 1)
        end_idx <- floor((i * total_pairs / total_cores))
        pairs <- sample_pairs[begin_idx:end_idx, ]
        out <- foreach(
              pair = iter(pairs, by = "row"), 
              .combine = rbind, 
              .verbose = verbose
            ) %do% {
            sample_a <- pair$sample_a
            sample_b <- pair$sample_b
            rxy <- ibdPair(
              list(dsmp[[sample_a]], dsmp[[sample_b]]), 
              c(coi[sample_a], coi[sample_b]), 
              afreq,
              M = 1, 
              pval = pval, 
              confreg = confint,
              rnull = rnull, 
              side = side, 
              alpha = alpha, 
              mnewton = mnewton,
              freqlog = TRUE, 
              reval = reval, 
              tol = tol, 
              logr = logr,
              neval = neval, 
              inull = inull, 
              nloc = nloc
            )
            estimate <- rxy$rhat
            p_value <- rxy$pval
            CI_lower <- range(rxy$confreg)[1]
            CI_upper <- range(rxy$confreg)[2]
            tibble::tibble(
              sample_a = sample_a, 
              sample_b = sample_b, 
              estimate = estimate, 
              p_value = p_value, 
              CI_lower = CI_lower, 
              CI_upper = CI_upper
            )

        }
        out
    }
    stopCluster(cl)
    setDefaultCluster(NULL)

    return(res)
}

#' Run `ibdEstM` in parallel
#'
#' This function uses `dcifer::ibdEstM` to estimate relatedness among 
#' all sample pairs in `dsmp`. `ibdEstM` differs from `dcifer::ibdPair` 
#' in that it will estimate the number of related strains (M) and 
#' return one estimate for each strain pair it estimates to be present.
#'
#' @inheritParams run_dcifer_bypop
#' @inheritParams dcifer::ibdDat
#' @param sample_pairs Tibble with two character columns, sample_a and 
#'   sample_b, specifying the sample combinations for which relatedness 
#'   should be estimated. If not provided, all possible sample pairs 
#'   will be run.
#' @param total_cores Integer specifying the number of cores to use.
#' @param verbose If TRUE, output from each parallel process will be 
#'   printed.
#'
#' @return A tibble with a character sample_a column, a character 
#'   sample_b column, an integer strain_pair column, and a double 
#'   estimate column. If `pval = TRUE`, there will also be a double 
#'   p_value column. If `confint = TRUE`, there will also be double 
#'   CI_lower and CI_upper columns.
#'
#' @details
#' If `rnull` is 0 or 1, Dcifer performs a one-sided hypothesis test, 
#' but it does a two-sided test if it is between 0 and 1. Therefore, if 
#' the scientific hypothesis of interest is, e.g., *r* > 0.25, it is 
#' necessary to divide the *p*-values returned by this function by two, 
#' and set *p*-values for *r* estimates below `rnull` to some 
#' arbitrarily high number if using something like the Benjamini-
#' Hochberg correction for multiple testing. These corrections are left 
#' to the user - they are NOT done in this function.
run_ibdestm <- function(
                        dsmp, 
                        coi, 
                        afreq, 
                        sample_pairs = NULL, 
                        pval = TRUE, 
                        confint = FALSE, 
                        rnull = 0, 
                        alpha = 0.05, 
                        total_cores = NULL, 
                        verbose = FALSE) {

    # Assemble ibdEstM inputs ------------------------------------------
    # Parameters used for all pairs
    nrs <- c(1e3, 1e2, 32, 16, 12, 10)
    revals <- mapply(generateReval, 1:6, nr = nrs)
    afreq <- lapply(afreq, log)
    nloc <- length(afreq)
    # Combine samples to get full set of pairs
    if (is.null(sample_pairs)) {
      sample_pairs_matrix <- combn(names(dsmp), 2)
      sample_pairs <- tibble(
        sample_a = sample_pairs_matrix[1,], 
        sample_b = sample_pairs_matrix[2,]
      )
    }

    # Determine side parameter based on rnull
    if (identical(rnull, 0L)) {
      side <- "right"
    } else if (identical(rnull, 1L)) {
      side <- "left"
    } else {
      side <- "two-sided"
    }

    # Set up cluster for parallel computing ----------------------------
    if (is.null(total_cores)) {
        total_cores <- parallelly::availableCores() - 1
    }
    if (is.null(getDefaultCluster())) {
        if (verbose) {
          cl <- makeCluster(total_cores, outfile = "")
        } else {
          cl <- makeCluster(total_cores)
        }
        setDefaultCluster(cl)
        registerDoParallel(cl)
    } else {
        cl <- getDefaultCluster()
    }

    # Run ibdEstM ------------------------------------------------------
    res <- foreach(
            i = 1:total_cores, 
            .combine = rbind, 
            .packages = c("dcifer", "foreach", "iterators")
          ) %dopar% {
        total_pairs <- nrow(sample_pairs)
        # Use floor to avoid off-by-one error when these equations 
        # don't yield a whole number. All pairs will still be included 
        # because end_idx will be a whole number when i = total_cores.
        begin_idx <- floor(((i - 1) * total_pairs / total_cores) + 1)
        end_idx <- floor((i * total_pairs / total_cores))
        pairs <- sample_pairs[begin_idx:end_idx, ]
        out <- foreach(
              pair = iter(pairs, by = "row"), 
              .combine = rbind, 
              .verbose = verbose
            ) %do% {
            sample_a <- pair$sample_a
            sample_b <- pair$sample_b
            rxy <- ibdEstM(
              list(dsmp[[sample_a]], dsmp[[sample_b]]), 
              c(coi[sample_a], coi[sample_b]), 
              afreq,
              Mmax = 6, 
              pval = pval, 
              confreg = confint,
              rnull = rnull, 
              side = side, 
              alpha = alpha, 
              equalr = FALSE, 
              freqlog = TRUE, 
              nrs = nrs, 
              revals = revals, 
              nloc = nloc
            )
            estimate <- rxy$rhat
            strain_pair <- 1:length(estimate)
            p_value <- rxy$pval
            CI_lower <- range(rxy$confreg)[1]
            CI_upper <- range(rxy$confreg)[2]
            tibble::tibble(
              sample_a = sample_a, 
              sample_b = sample_b, 
              estimate = estimate, 
              strain_pair = strain_pair, 
              p_value = p_value, 
              CI_lower = CI_lower, 
              CI_upper = CI_upper
            )

        }
        out
    }
    stopCluster(cl)
    setDefaultCluster(NULL)

    return(res)
}

#' Reformat Dcifer output and write to disk
#'
#' This function receives a tibble containing Dcifer results, lightly 
#' reformats it, and saves it to the provided path.
#'
#' @param dcifer_results A tibble containing character sample_a and 
#'   sample_b columns, a double estimate column, and optionally double 
#'   p_value, CI_lower, and CI_upper columns. This is the output of 
#'   `run_ibdpair`.
#' @param out_path Path to save Dcifer results as a TSV.
write_dcifer_output <- function(dcifer_results, out_path) {
  dcifer_results %>%
    rename(
      specimen_name_a = sample_a, 
      specimen_name_b = sample_b, 
      btwn_host_rel = estimate
    ) %>%
    write_tsv(out_path)
}

set.seed(arg$seed)

# Read data and/or calculate COI and allele frequencies ----------------
# Read in alleles to Dcifer format
allele_table <- create_allele_table_input(
  arg$allele_table, 
  specimen_name_col = arg$specimen_name_col, 
  target_name_col = arg$target_name_col, 
  target_value_col = arg$target_value_col
)
dcifer_alleles <- dcifer::formatDat(
  allele_table, 
  svar = "specimen_name", 
  lvar = "target_name", 
  avar = "target_value"
)
# If no COI input was provided, use Dcifer's built-in naive estimation
if (is.null(arg$coi_table)) {
  coi <- dcifer::getCOI(dcifer_alleles)
  # Remove locus used for estimation from vector names
  names(coi) <- names(coi) %>%
    str_split_i("\\.", 1)
} else {
  coi <- create_coi_input(
    arg$coi_table, 
    dcifer_alleles, 
    specimen_name_col = arg$specimen_name_col
  )
}
# If no allele frequencies were provided, use Dcifer's built-in naive 
# estimation
if (is.null(arg$allele_freq_table)) {
  # Compute for all data at once if pop_name_col not provided
  if (is.null(arg$pop_name_col)) {
    allele_freqs <- dcifer::calcAfreq(dcifer_alleles, coi, tol = 1e-5)
  # If pop_name_col was included, check for and read in specimen_metadata
  } else {
    if (is.null(arg$specimen_metadata)) {
      stop("pop_name_col provided but specimen_metadata missing", call. = FALSE)
    }
    # Read in specimen metadata table
    specimen_metadata <- create_specimen_metadata_input(
        arg$specimen_metadata, 
        coi, 
        specimen_name_col = arg$specimen_name_col, 
        pop_name_col = arg$pop_name_col
      )
  }
# If allele frequencies were provided, read them in
} else {
  allele_freqs <- create_allele_freq_input(
    arg$allele_freq_table, 
    dcifer_alleles, 
    target_name_col = arg$target_name_col, 
    target_value_col = arg$target_value_col
  )
}

# Run Dcifer and save output -------------------------------------------
if (is.null(arg$pop_name_col)) {
  if (arg$use_estm) {
    dcifer_res <- run_ibdestm(
        dcifer_alleles, 
        coi, 
        allele_freqs, 
        confint = TRUE, 
        rnull = arg$rnull, 
        alpha = arg$alpha, 
        total_cores = arg$threads, 
        verbose = arg$verbose
      )
  } else {
    dcifer_res <- run_ibdpair(
        dcifer_alleles, 
        coi, 
        allele_freqs, 
        confint = TRUE, 
        rnull = arg$rnull, 
        alpha = arg$alpha, 
        total_cores = arg$threads, 
        verbose = arg$verbose
      )
  }
} else {
  dcifer_res <- run_dcifer_bypop(
    arg$allele_table, 
    coi, 
    specimen_metadata, 
    specimen_name_col = arg$specimen_name_col, 
    target_name_col = arg$target_name_col, 
    target_value_col = arg$target_value_col, 
    use_estm = arg$use_estm, 
    confint = TRUE, 
    rnull = arg$rnull, 
    alpha = arg$alpha, 
    total_cores = arg$threads, 
    verbose = arg$verbose
  )
}
write_dcifer_output(dcifer_res, arg$btwn_host_rel_output)
