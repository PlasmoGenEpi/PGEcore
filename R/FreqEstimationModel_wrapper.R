#' Average COI from a table path or a numeric value
#'
#' @keywords internal
calculate_avg_COI <- function(coi_path) {
  coi <- suppressWarnings(as.numeric(coi_path))
  if (!is.na(coi)) {
    return(as.numeric(coi_path))
  }
  COI_table <- readr::read_tsv(
    coi_path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      coi = readr::col_number()
    )
  )
  rules <- validate::validator(
    !is.na(specimen_name),
    !is.na(coi)
  )
  stop_on_validate_fails(COI_table, rules, "coi")
  mean(COI_table$coi)
}

#' Read FEM loci-group definitions
#'
#' @keywords internal
read_fem_groups <- function(groups_path) {
  readr::read_tsv(
    groups_path,
    col_types = readr::cols(
      group_id = readr::col_character(),
      gene_id = readr::col_character(),
      aa_position = readr::col_integer()
    )
  )
}

#' Keep bi-allelic (and report mono-allelic) amino-acid targets
#'
#' @keywords internal
check_fem_biallelic <- function(input_data) {
  mutants <- input_data |>
    dplyr::filter(.data$ref_aa != .data$aa) |>
    dplyr::distinct(.data$unique_targets, .data$aa, .keep_all = TRUE)
  bad_targets <- mutants |>
    dplyr::filter(.data$unique_targets %in% duplicated(.data$unique_targets))
  if (nrow(bad_targets) > 0) {
    warning("Not biallelic, dropped", bad_targets, call. = FALSE)
  }
  monoallelic <- input_data |>
    dplyr::group_by(.data$unique_targets) |>
    dplyr::filter(dplyr::n_distinct(.data$aa) == 1) |>
    dplyr::distinct(
      .data$gene_id, .data$aa_position, .data$aa, .data$unique_targets
    ) |>
    dplyr::ungroup()

  only_biallelic <- input_data |>
    dplyr::filter(
      !(.data$unique_targets %in% bad_targets$unique_targets) &
        !(.data$unique_targets %in% monoallelic$unique_targets)
    )

  alt_calls <- dplyr::filter(only_biallelic, .data$ref_aa != .data$aa)
  list(only_biallelic, alt_calls, monoallelic)
}

#' Read amino-acid calls for FEM
#'
#' @keywords internal
read_fem_aa_calls <- function(input_path) {
  input_data <- readr::read_tsv(
    input_path,
    col_types = readr::cols(
      specimen_name = readr::col_character(),
      aa_position = readr::col_integer(),
      target_name = readr::col_character(),
      gene_id = readr::col_character(),
      ref_aa = readr::col_character(),
      aa = readr::col_character()
    )
  )
  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(gene_id),
    is.integer(aa_position),
    is.character(ref_aa),
    is.character(aa),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(ref_aa),
    !is.na(aa)
  )
  stop_on_validate_fails(input_data, rules, "aa_calls")
  input_data
}

#' Build the genotype matrix and related FEM inputs for one group
#'
#' @keywords internal
create_FEM_input <- function(input_data, groups, group_id) {
  input_data$unique_targets <- paste(
    input_data$gene_id, input_data$aa_position, sep = ":"
  )
  groups <- groups[groups$group_id == group_id, ]
  group_targets <- paste(groups$gene_id, groups$aa_position, sep = ":")
  input_data <- input_data[input_data$unique_targets %in% group_targets, ]
  data_list <- check_fem_biallelic(input_data)
  input_data <- data_list[[1]]
  alt_alleles <- data_list[[2]]
  monos <- data_list[[3]]

  unique_targets <- unique(input_data$unique_targets)
  unique_sample_ids <- unique(input_data$specimen_name)

  sample_matrix <- matrix(
    99,
    nrow = length(unique_sample_ids),
    ncol = length(unique_targets)
  )
  colnames(sample_matrix) <- unique_targets
  rownames(sample_matrix) <- unique_sample_ids

  for (sample in input_data$specimen_name) {
    cut_df <- input_data[input_data$specimen_name == sample, ]
    for (unique_target in unique_targets) {
      cut_cut_df <- cut_df[cut_df$unique_targets == unique_target, ]
      nvals <- length(unique(cut_cut_df$aa))
      if (nvals > 2) {
        stop("Too many alleles for FEM", call. = FALSE)
      }
      if (nvals == 2) {
        sample_matrix[sample, unique_target] <- 0.5
      }
      if (nvals == 1) {
        if (cut_cut_df[1, "ref_aa"] == cut_cut_df[1, "aa"]) {
          val <- 0
        } else {
          val <- 1
        }
        sample_matrix[sample, unique_target] <- val
      }
    }
  }
  list(sample_matrix, alt_alleles, nrow(sample_matrix), monos)
}

#' Run FreqEstimationModel MCMC for one group
#'
#' @param sample_matrix_list Output of [create_FEM_input()].
#' @param COI Average complexity of infection.
#' @param threads Number of threads.
#' @param seed Random seed.
#' @param num_chains Number of MCMC chains to run. At least two are needed to
#'   compute the Gelman-Rubin R-hat convergence diagnostic.
#'
#' @return A list with the population frequency table, MCMC runtime, marker
#'   names, alternate alleles, group size, mono-allelic loci, and a
#'   `convergence_diag` data frame.
#' @keywords internal
run_FreqEstimationModel <- function(sample_matrix_list,
                                    COI,
                                    threads,
                                    seed,
                                    num_chains = 3L) {
  check_suggested_pkgs(
    c("FreqEstimationModel", "plyr", "coda", "abind", "foreach", "doMC"),
    "FreqEstimationModel MCMC"
  )
  check_suggested_pkg("posterior", "FEM convergence diagnostics")

  sample_matrix <- sample_matrix_list[[1]]
  alt_alleles <- sample_matrix_list[[2]]
  num_group <- sample_matrix_list[[3]]
  monos <- sample_matrix_list[[4]]
  data_summary <- list()
  data_summary$Data <- sample_matrix
  runtime <- system.time({
    thinning_interval <- 1
    no_traces_preburnin <- 10000
    no_mcmc_chains <- num_chains
    NGS <- FALSE
    log_like_zero <- FALSE
    mcmc_variable_list <- list(
      no_mcmc_chains = no_mcmc_chains,
      no_traces_preburnin = no_traces_preburnin,
      thinning_interval = thinning_interval,
      NGS = NGS,
      log_like_zero = log_like_zero
    )
    if (!NGS | log_like_zero) {
      augment_missing_data <- TRUE
    } else {
      augment_missing_data <- FALSE
    }
    if (log_like_zero) {
      moi_prior <- "Uniform"
    } else {
      moi_prior <- "Poisson"
    }
    moi_max <- 8
    moi_hyperparameter <- COI
    moi_size_hyperparameter <- 0.5
    moi_prior_min2 <- NULL
    moi_initial <- NULL
    moi_list <- list(
      moi_hyperparameter = moi_hyperparameter,
      moi_size_hyperparameter = moi_size_hyperparameter,
      moi_prior = moi_prior,
      moi_max = moi_max,
      moi_prior_min2 = moi_prior_min2,
      moi_initial = moi_initial
    )
    processed_data_list <- FreqEstimationModel::preprocess_data(
      data_summary,
      log_like_zero,
      NGS,
      augment_missing_data,
      moi_prior_min2
    )
    frequency_hyperparameter <- rep(1, processed_data_list$no_haplotypes)
    frequency_initial <- NULL
    frequency_list <- list(
      frequency_hyperparameter = frequency_hyperparameter,
      frequency_initial = frequency_initial
    )

    set.seed(as.numeric(seed))
    results <- FreqEstimationModel::mcmc_sampling_parallel(
      processed_data_list,
      moi_list,
      frequency_list,
      mcmc_variable_list,
      cores_max = threads
    )

    burnin <- 1:(0.5 * mcmc_variable_list$no_traces_preburnin)
    if (mcmc_variable_list$no_mcmc_chains > 1) {
      alply_genotype_freq_store_chains_burnin <- plyr::alply(
        results$genotype_freq_store_chains[-burnin, , ],
        3
      )
    } else {
      alply_genotype_freq_store_chains_burnin <-
        results$genotype_freq_store_chains[-burnin, , ]
    }

    mcmc_frequency_chains <- coda::mcmc.list(lapply(
      alply_genotype_freq_store_chains_burnin,
      coda::mcmc,
      start = (max(burnin) + 1) * mcmc_variable_list$thinning_interval,
      end = mcmc_variable_list$no_traces_preburnin *
        mcmc_variable_list$thinning_interval,
      thin = mcmc_variable_list$thinning_interval
    ))

    mcmc_As <- abind::abind(
      plyr::alply(results$genotype_count_store_chains[-burnin, , , ], 4),
      along = 1
    )
    mcmc_mois <- apply(mcmc_As, c(1, 2), sum)

    pop_freq <- cbind(
      freq = summary(mcmc_frequency_chains)$statistics[, "Mean"],
      median_freq = summary(mcmc_frequency_chains)$quantiles[, 3],
      "CI_2.5" = summary(mcmc_frequency_chains)$quantiles[, 1],
      "CI_97.5" = summary(mcmc_frequency_chains)$quantiles[, 5]
    )
    pop_prev <- 1 - (1 - pop_freq)^median(mcmc_mois)
    inf_prev <- t(apply(mcmc_As, 2, function(x) colMeans(x > 0)))
    prev <- colMeans(inf_prev)
    pop_freq <- cbind(pop_freq, prev)
  })
  sequence_column <- data.frame(
    sequence = rownames(pop_freq),
    stringsAsFactors = FALSE
  )
  pop_freq <- cbind(sequence_column, pop_freq)
  rownames(pop_freq) <- NULL

  # One variable per haplotype frequency chain
  convergence_diag <- summarize_convergence_draws(
    posterior::as_draws_array(mcmc_frequency_chains)
  )

  list(
    plsf_table = pop_freq,
    runtime = runtime,
    names = processed_data_list[["markerID"]],
    alt_allele = alt_alleles,
    num_group = num_group,
    monos = monos,
    convergence_diag = convergence_diag
  )
}

#' Convert a FEM binary haplotype string to a STAVE variant
#'
#' @keywords internal
bin2STAVE <- function(chars, names, alt_alleles, monos) {
  check_suggested_pkg("variantstring", "STAVE strings via bin2STAVE()")
  char_ix <- 1
  chars_split <- strsplit(chars, "")[[1]]
  long_form <- tibble::tibble(
    gene = character(),
    pos = numeric(),
    n_aa = numeric(),
    het = logical(),
    phased = logical(),
    aa = character(),
    read_count = numeric()
  )

  for (char in chars_split) {
    if (char == 1) {
      call <- "aa"
    } else {
      call <- "ref_aa"
    }
    alt_current <- alt_alleles[alt_alleles$unique_targets == names[char_ix], ]
    gene_current <- stringr::str_split(names[char_ix], ":")[[1]][1]
    pos_current <- stringr::str_split(names[char_ix], ":")[[1]][2]
    alt_current <- alt_current[1, call]
    alt_current <- as.character(alt_current)
    char_ix <- char_ix + 1
    long_form <- tibble::add_row(
      long_form,
      gene = gene_current,
      pos = as.numeric(pos_current),
      n_aa = NA,
      het = NA,
      phased = TRUE,
      aa = alt_current,
      read_count = NA
    )
  }
  mono_df <- monos |>
    dplyr::rename(gene = "gene_id", pos = "aa_position") |>
    dplyr::mutate(
      n_aa = NA,
      het = FALSE,
      phased = TRUE,
      read_count = NA
    ) |>
    dplyr::select(-"unique_targets")
  long_form <- rbind(long_form, mono_df)
  long_form <- long_form |>
    dplyr::group_by(.data$gene, .data$pos) |>
    dplyr::mutate(n_aa = dplyr::n())
  long_form <- long_form |>
    dplyr::mutate(het = ifelse(.data$n_aa > 1, TRUE, FALSE))
  variantstring::long_to_variant(list(long_form))
}

#' Format FEM output for one group
#'
#' @keywords internal
format_single_group_output <- function(pop_freq_list) {
  input_list <- pop_freq_list[[1]]
  names <- pop_freq_list[[3]]
  alt_alleles <- pop_freq_list[[4]]
  num_group <- pop_freq_list[[5]]
  monos <- pop_freq_list[[6]]

  input_list |>
    dplyr::mutate(
      variant = purrr::map_chr(
        .data$sequence,
        bin2STAVE,
        names,
        alt_alleles,
        monos
      ),
      sample_total = num_group
    ) |>
    tibble::as_tibble()
}

#' Format output for an invariant FEM group
#'
#' @keywords internal
format_invariant_group_output <- function(aa_calls, groups, group) {
  check_suggested_pkg(
    "variantstring",
    "STAVE strings via format_invariant_group_output()"
  )
  group_loci <- groups |>
    dplyr::filter(.data$group_id == group) |>
    tidyr::unite("loci_name", "gene_id", "aa_position", sep = ":") |>
    dplyr::pull("loci_name")
  group_calls <- aa_calls |>
    dplyr::distinct(
      .data$specimen_name, .data$gene_id, .data$aa_position, .data$aa
    ) |>
    tidyr::unite(
      "loci_name", "gene_id", "aa_position",
      sep = ":",
      remove = FALSE
    ) |>
    dplyr::filter(.data$loci_name %in% group_loci)
  num_group <- dplyr::n_distinct(group_calls$specimen_name)
  n_alleles_per_locus <- group_calls |>
    dplyr::group_by(.data$gene_id, .data$aa_position) |>
    dplyr::mutate(n_allele = dplyr::n_distinct(.data$aa)) |>
    dplyr::ungroup()
  if (any(n_alleles_per_locus$n_allele > 1)) {
    stop(
      "Group ",
      group,
      " is presumed invariant but some loci have multiple alleles. FEM ",
      "inputs are empty for some other reason.",
      call. = FALSE
    )
  }
  variant_stave <- group_calls |>
    dplyr::distinct(.data$gene_id, .data$aa_position, .data$aa) |>
    dplyr::rename(gene = "gene_id", pos = "aa_position") |>
    dplyr::group_by(.data$gene, .data$pos) |>
    dplyr::mutate(n_aa = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      het = FALSE,
      phased = TRUE,
      read_count = NA
    ) |>
    dplyr::relocate("aa", .before = "read_count") |>
    list() |>
    variantstring::long_to_variant()

  tibble::tibble(
    sequence = NA,
    freq = 1,
    median_freq = 1,
    CI_2.5 = 1,
    CI_97.5 = 1,
    prev = 1,
    variant = variant_stave,
    sample_total = num_group
  )
}

#' Estimate multilocus allele frequencies with FreqEstimationModel
#'
#' File-oriented entry point used by the `FreqEstimationModel_wrapper` CLI.
#' Optional **FreqEstimationModel**, **variantstring**, **posterior**, and
#' parallel helpers (**foreach**, **doMC**, plus **plyr**, **coda**, **abind**)
#' must be installed separately.
#'
#' @param aa_calls Path to amino-acid call TSV.
#' @param coi Path to COI TSV, or a numeric average COI.
#' @param groups Path to group TSV (`group_id`, `gene_id`, `aa_position`).
#' @param mlaf_output Output TSV path.
#' @param threads Number of threads.
#' @param seed Random seed.
#' @param num_chains Number of MCMC chains to run per group. At least two are
#'   needed to compute the Gelman-Rubin R-hat convergence diagnostic.
#' @param convergence_output Output TSV path for per-group MCMC convergence
#'   diagnostics, with the columns `group_id`, `variable`, `mean`, `median`,
#'   `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`.
#'
#' @return The formatted output data frame (also written to `mlaf_output`).
#' @export
FreqEstimationModel_wrapper <- function(aa_calls,
                                        coi,
                                        groups,
                                        mlaf_output,
                                        threads = 1L,
                                        seed = 1L,
                                        num_chains = 3L,
                                        convergence_output = "convergence_diag.tsv") {
  check_suggested_pkg(
    "FreqEstimationModel",
    "multilocus frequencies via FreqEstimationModel_wrapper()"
  )
  check_suggested_pkg(
    "variantstring",
    "STAVE strings via FreqEstimationModel_wrapper()"
  )
  check_suggested_pkg(
    "posterior",
    "MCMC convergence diagnostics via FreqEstimationModel_wrapper()"
  )
  check_suggested_pkgs(
    c("foreach", "doMC", "plyr", "coda", "abind"),
    "FreqEstimationModel MCMC helpers"
  )

  if (is.null(aa_calls) || is.null(coi) || is.null(groups) || is.null(mlaf_output)) {
    stop(
      "--aa_calls, --coi, --groups, and --mlaf_output are required",
      call. = FALSE
    )
  }

  aa_tbl <- read_fem_aa_calls(aa_calls)
  groups_tbl <- read_fem_groups(groups)
  COI <- calculate_avg_COI(coi)
  overall_output <- data.frame(
    sequence = character(),
    freq = numeric(),
    median_freq = numeric(),
    CI_2.5 = numeric(),
    CI_97.5 = numeric()
  )
  # Invariant groups skip MCMC and contribute no diagnostic rows.
  convergence_diags <- list()

  for (group in unique(groups_tbl$group_id)) {
    fem_input <- create_FEM_input(aa_tbl, groups_tbl, group)
    if (sum(dim(fem_input[[1]])) == 0) {
      fem_plsf <- format_invariant_group_output(aa_tbl, groups_tbl, group)
    } else {
      fem_results <- run_FreqEstimationModel(
        fem_input,
        COI,
        threads,
        seed,
        num_chains
      )
      fem_plsf <- format_single_group_output(fem_results)
      convergence_diags[[group]] <- fem_results$convergence_diag |>
        dplyr::mutate(group_id = group) |>
        dplyr::relocate("group_id")
    }
    fem_plsf$group_id <- group
    overall_output <- rbind(overall_output, fem_plsf)
  }

  overall_output <- apply(overall_output, 2, as.character)
  overall_output_df <- data.frame(overall_output)
  readr::write_tsv(overall_output_df, mlaf_output)

  readr::write_tsv(dplyr::bind_rows(convergence_diags), convergence_output)

  overall_output_df
}
