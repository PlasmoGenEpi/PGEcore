# Estimate IBD-based relatedness with Dcifer

# Load required libraries ----------------------------------------------
library(dcifer)
# These will be referenced without the `package::` construct, and thus 
# are loaded second to avoid masking
library(doParallel)
library(dplyr)
library(magrittr)
library(optparse)
library(parallel)
library(parallelly)
library(purrr)
library(readr)
library(stringr)
library(tibble)
library(tidyr)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option(
    "--allele_table", 
    help = str_c(
      "TSV containing alleles, with the columns: specimen_id, target_id, ", 
      "read_count, and seq. Required."
    )
  ), 
  make_option(
    "--coi_table", 
    help = 
      str_c(
        "TSV containing specimen COIs, with the columns: specimen_id and ", 
        "coi. Optional."
      )
  ), 
  make_option(
    "--allele_freq_table", 
    help = 
      str_c(
        "TSV containing single locus allele frequencies, with the columns: ", 
        "gene_id, seq, freq, and total. Optional."
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
      "specimen_id, btwn_host_rel. Required."
    )
  )
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
if (interactive()) {
  arg$allele_table <- "../../data/example_allele_table.tsv"
  arg$coi_table <- "../../data/example_coi_table.tsv"
  arg$btwn_host_rel_output <- "../../btwn_host_rel.tsv"
}

#' Read allele table into a tibble
#'
#' Read the allele table TSV into a tibble and then call 
#' `dcifer::formatDat` to convert into list format.
#'
#' @param allele_table_path Path to the TSV file containing the allele 
#'   table. There should be character columns for specimen_id, 
#'   target_id, and seq.
#'
#' @return The list format output by `dcifer::readDat` and 
#'   `dcifer::formatDat`.
create_allele_table_input <- function(allele_table_path) {

  # Read in table
  allele_table <- read_tsv(
    allele_table_path, 
    col_types = cols(.default = col_character()), 
    progress = FALSE
  )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_id), 
    is.character(target_id), 
    is.character(seq), 
    ! is.na(specimen_id), 
    ! is.na(target_id), 
    ! is.na(seq)
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

  # Convert to Dcifer format and return
  dcifer::formatDat(
    allele_table, 
    svar = "specimen_id", 
    lvar = "target_id", 
    avar = "seq"
  )

}

#' Read COI table into format needed by Dcifer
#'
#' Read the TSV of specimen COIs into a tibble, join with specimen IDs 
#' from the allele list to ensure order matches, and return as a vector.
#'
#' @param coi_path Path to TSV containing specimen COIs. It should 
#'   have a character specimen_id column and an integer coi column.
#' @param allele_list The list format output by `dcifer::readDat` and 
#'   `dcifer::formatDat`.
#'
#' @return Vector of COI values, one for each sample.
create_coi_input <- function(coi_path, allele_list) {

  # Read input table
  coi <- read_tsv(
    coi_path, 
    col_types = cols(.default = col_character(), coi = col_integer()), 
    progress = FALSE
  )

  # Validate fields
  rules <- validate::validator(
    is.character(specimen_id), 
    is.integer(coi), 
    ! is.na(specimen_id), 
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
  specimen_inallele_notincoi <- setdiff(names(allele_list), coi$specimen_id)
  if (length(specimen_inallele_notincoi) > 0) {
    stop(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_inallele_notincoi, collapse = " "), 
      call. = FALSE
    )
  }
  specimen_incoi_notinallele <- setdiff(coi$specimen_id, names(allele_list))
  if (length(specimen_incoi_notinallele) > 0) {
    warning(
      "The following specimen IDs appear in the allele table and not in the ", 
      "COI table: ", 
      str_c(specimen_incoi_notinallele, collapse = " "), 
      call. = FALSE
    )
  }

  # Join COI to allele list to ensure order matches
  coi <- tibble(specimen_id = names(allele_list)) %>%
    left_join(coi, by = "specimen_id")

  # Create and return named vector of COI
  setNames(coi$coi, coi$specimen_id)

}

#' Read allele frequency into a list format
#'
#' Read a TSV of allele frequencies into a tibble, then reformat into 
#' the list format taken by 'dcifer::ibdDat` and `dcifer::ibdPair`.
#'
#' @param allele_freq_path Path to the TSV file containing the allele 
#'   frequencies. There should be character columns for target_id and 
#'   seq, and a double freq column.
#' @inheritParams create_allele_table_input
#'
#' @return A named list with target_id as the names. The values are a 
#'   named vector with seq as the name and freq as the value. This is 
#'   format taken by `dcifer::ibdDat` and `dcifer::ibdPair`.
create_allele_freq_input <- function(allele_freq_path, allele_list) {

  # Read in table
  allele_freqs <- read_tsv(
    allele_freq_path, 
    col_types = cols(
      .default = col_character(), 
      freq = col_double()
    ), 
    progress = FALSE
  )

  # Validate fields
  rules <- validate::validator(
    is.character(target_id), 
    is.character(seq), 
    is.double(freq), 
    ! is.na(target_id), 
    ! is.na(seq), 
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
    unite(allele, target_id, seq, sep = ":") %$%
    allele
  allele_table_alleles <- list_c(allele_list) %>%
    tibble(target_id = names(.), seqs = .) %>%
    mutate(seqs = map(seqs, names)) %>%
    unnest(seqs) %>%
    unite(alleles, target_id, seqs, sep = ":") %$%
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
    lvar = "target_id", 
    avar = "seq", 
    fvar = "freq"
  )

}

#' Run Dcifer in parallel
#'
#' This function uses `dcifer::ibdPair` to estimate relatedness among 
#' all sample pairs in `dsmp`. It essentially replicates the 
#' functionality of `dcifer::ibdDat` while allowing for parallel 
#' execution.
#'
#' @inheritParams dcifer::ibdDat
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
#' (very) lightly modified by Alfred Hubbard.
run_dcifer <- function(
                            dsmp, 
                            coi, 
                            afreq, 
                            dsmp2 = NULL, 
                            coi2 = NULL, 
                            pval = TRUE, 
                            confint = FALSE, 
                            rnull = 0, 
                            alpha = 0.05, 
                            nr = 1000, 
                            reval = NULL, 
                            total_cores = NULL, 
                            verbose = FALSE) {

    dwithin <- is.null(dsmp2)
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
    nsmp <- length(dsmp)
    snames <- names(dsmp)
    if (dwithin) {
        dsmp2 <- dsmp
        coi2 <- coi
    }
    nsmp2 <- length(dsmp2)
    snames2 <- names(dsmp2)

    sample_pairs <- expand.grid(1:nsmp, 1:nsmp2) |>
        dplyr::filter(Var1 < Var2)

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
            ix <- pair$Var1
            iy <- pair$Var2
            rxy <- ibdPair(list(dsmp[[ix]], dsmp2[[iy]]), c(
                coi[ix],
                coi2[iy]
            ), afreq,
            M = 1, pval = pval, confreg = confint,
            rnull = rnull, alpha = alpha, mnewton = mnewton,
            freqlog = TRUE, reval = reval, tol = tol, logr = logr,
            neval = neval, inull = inull, nloc = nloc
            )
            estimate <- rxy$rhat
            p_value <- rxy$pval
            CI_lower <- range(rxy$confreg)[1]
            CI_upper <- range(rxy$confreg)[2]
            tibble::tibble(
              sample_a = names(dsmp)[ix], 
              sample_b = names(dsmp)[iy], 
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

#' Reformat Dcifer output and write to disk
#'
#' This function receives a tibble containing Dcifer results, lightly 
#' reformats it, and saves it to the provided path.
#'
#' @param dcifer_results A tibble containing character sample_a and 
#'   sample_b columns, a double estimate column, and optionally double 
#'   p_value, CI_lower, and CI_upper columns. This is the output of 
#'   `run_dcifer`.
#' @param out_path Path to save Dcifer results as a TSV.
write_dcifer_output <- function(dcifer_results, out_path) {
  dcifer_results %>%
    rename(
      specimen_id_a = sample_a, 
      specimen_id_b = sample_b, 
      btwn_host_rel = estimate
    ) %>%
    write_tsv(out_path)
}

set.seed(arg$seed)

# Read data and/or calculate COI and allele frequencies ----------------
dcifer_alleles <- create_allele_table_input(arg$allele_table)
# If no COI input was provided, use Dcifer's built-in naive estimation
if (is.null(arg$coi_table)) {
  coi <- dcifer::getCOI(dcifer_alleles)
} else {
  coi <- create_coi_input(arg$coi_table, dcifer_alleles)
}
# If no allele frequencies were provided, use Dcifer's built-in naive 
# estimation
if (is.null(arg$allele_freq_table)) {
  allele_freqs <- dcifer::calcAfreq(dcifer_alleles, coi, tol = 1e-5)
} else {
  allele_freqs <- create_allele_freq_input(
    arg$allele_freq_table, 
    dcifer_alleles
  )
}

# Compute relatedness and save -----------------------------------------
run_dcifer(
    dcifer_alleles, 
    coi, 
    allele_freqs, 
    confint = TRUE, 
    rnull = arg$rnull, 
    alpha = arg$alpha, 
    total_cores = arg$threads, 
    verbose = arg$verbose
  ) %>%
  write_dcifer_output(arg$btwn_host_rel_output)
