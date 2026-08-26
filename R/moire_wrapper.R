#' Create a MOIRe input object from an allele table
#'
#' Reads a TSV of allele presence, validates columns, and packages MCMC
#' parameters for [run_moire()].
#'
#' @keywords internal
create_moire_input <- function(input_path,
                               allow_relatedness,
                               burnin,
                               samples_per_chain,
                               thin,
                               verbose,
                               eps_pos_alpha,
                               eps_pos_beta,
                               eps_neg_alpha,
                               eps_neg_beta,
                               r_alpha,
                               r_beta,
                               mean_coi_shape,
                               mean_coi_scale,
                               max_eps_pos,
                               max_eps_neg,
                               record_latent_genotypes,
                               pt_chains,
                               pt_grad_lower,
                               pt_num_threads,
                               adapt_temp,
                               max_runtime,
                               n_chains,
                               threads) {
  check_suggested_pkg("checkmate", "MOIRe input validation")

  message("Reading input data")
  input_data <- utils::read.csv(
    input_path,
    na.strings = "NA",
    sep = "\t",
    colClasses = c(specimen_name = "character")
  )

  validate_required_columns(
    input_data,
    c("specimen_name", "target_name", "seq"),
    input_path
  )

  message("Validating input format")
  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(seq),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(seq)
  )
  message("Confronting input data with validation rules")
  stop_on_validate_fails(input_data, rules, input_path)

  moire_data <- input_data |>
    dplyr::select("specimen_name", "target_name", "seq") |>
    dplyr::rename(
      sample_id = "specimen_name",
      locus = "target_name",
      allele = "seq"
    )

  message("Creating Moire object")

  # This wrapper does not support passing MOIRe's pt_grad argument. Instead a
  # uniformly spaced sequence of rungs down to pt_grad_lower is generated,
  # because the highly tempered distributions often do not swap well. The
  # ladder must be descending (cold chain, temperature 1.0, first): MOIRe
  # reports samples from rung 1 as the cold chain, so an ascending ladder makes
  # it sample the pure-prior rung and corrupts the posterior summaries.
  if (pt_chains > 1) {
    pt_chains <- seq(from = 1, to = pt_grad_lower, length.out = pt_chains)
  } else {
    pt_chains <- 1
  }

  moire_object <- list(
    moire_data = moire_data,
    moire_parameters = list(
      allow_relatedness = allow_relatedness,
      burnin = burnin,
      samples_per_chain = samples_per_chain,
      thin = thin,
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
      max_runtime = max_runtime,
      num_chains = n_chains,
      num_cores = threads
    )
  )

  p <- moire_object$moire_parameters
  checkmate::assert_logical(p$allow_relatedness, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$burnin, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$samples_per_chain, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$thin, any.missing = FALSE, len = 1)
  checkmate::assert_logical(p$verbose, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$eps_pos_alpha, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$eps_pos_beta, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$eps_neg_alpha, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$eps_neg_beta, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$r_alpha, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$r_beta, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$mean_coi_shape, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$mean_coi_scale, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$max_eps_pos, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$max_eps_neg, any.missing = FALSE, len = 1)
  checkmate::assert_logical(p$record_latent_genotypes, any.missing = FALSE, len = 1)
  # pt_grad_lower only shapes the pt_chains ladder and is not stored in
  # moire_parameters, so the incoming argument is validated directly.
  checkmate::assert_numeric(pt_grad_lower, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$pt_num_threads, any.missing = FALSE, len = 1)
  checkmate::assert_logical(p$adapt_temp, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$max_runtime, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$num_chains, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$num_cores, any.missing = FALSE, len = 1)

  message("Returning Moire object")
  moire_object
}

#' Run MOIRe MCMC analysis
#'
#' Runs MOIRe MCMC on a prepared `moire_object`. Requires **moire** (Suggests).
#' For reading allele tables and writing summary TSVs, use [moire_wrapper()].
#'
#' ## Inputs
#'
#' - **`moire_object`**: List with `moire_data` and `moire_parameters`, as
#'   created by [create_moire_input()] (via [moire_wrapper()]).
#'
#' ## Outputs
#'
#' - Returns the object from [moire::run_mcmc()] (not written to disk).
#'
#' ## Running
#'
#' ```r
#' run_moire(moire_object)
#' ```
#'
#' File and CLI users should call [moire_wrapper()] /
#' `Rscript exec/moire_wrapper ...`.
#'
#' Requires **moire** (Suggests).
#'
#' @param moire_object List with `moire_data` and `moire_parameters`. See
#'   *Inputs*.
#'
#' @return The object returned by [moire::run_mcmc()].
#'
#' @seealso [moire_wrapper()], `vignette("input-formats", package = "PGEcore")`
#'
#' @export
run_moire <- function(moire_object) {
  check_suggested_pkg("moire", "MCMC analysis via run_moire()")

  moire_data <- moire::load_long_form_data(moire_object$moire_data)
  moire_parameters <- moire_object$moire_parameters

  moire::run_mcmc(
    moire_data,
    moire_data$is_missing,
    allow_relatedness = moire_parameters$allow_relatedness,
    burnin = moire_parameters$burnin,
    samples_per_chain = moire_parameters$samples_per_chain,
    thin = moire_parameters$thin,
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
    max_runtime = moire_parameters$max_runtime,
    num_chains = moire_parameters$num_chains,
    num_cores = moire_parameters$num_cores
  )
}

#' Summarize MOIRe MCMC results and write TSV files
#'
#' @keywords internal
summarize_and_write_moire_results <- function(moire_object,
                                              mcmc_results,
                                              coi_output,
                                              he_output,
                                              allele_freq_output,
                                              relatedness_output,
                                              effective_coi_output) {
  check_suggested_pkg("moire", "summarizing MOIRe MCMC results")
  check_suggested_pkg("checkmate", "MOIRe summary checks")

  coi_summary <- moire::summarize_coi(mcmc_results) |>
    dplyr::rename(specimen_name = "sample_id", coi = "post_coi_mean")
  he_summary <- moire::summarize_he(mcmc_results) |>
    dplyr::rename(target_name = "locus", he = "post_stat_mean")
  allele_freq_summary <- moire::summarize_allele_freqs(mcmc_results) |>
    dplyr::rename(
      target_name = "locus",
      freq = "post_allele_freqs_mean",
      seq = "allele"
    )
  relatedness_summary <- moire::summarize_relatedness(mcmc_results) |>
    dplyr::rename(
      specimen_name = "sample_id",
      within_host_rel = "post_relatedness_mean"
    )
  effective_coi_summary <- moire::summarize_effective_coi(mcmc_results) |>
    dplyr::rename(
      specimen_name = "sample_id",
      ecoi = "post_effective_coi_mean"
    )

  missing_target_names <- unique(moire_object$moire_data$locus[
    !moire_object$moire_data$locus %in% allele_freq_summary$target_name
  ])
  present_target_names <- unique(moire_object$moire_data$locus[
    moire_object$moire_data$locus %in% allele_freq_summary$target_name
  ])

  checkmate::assert(
    length(unique(moire_object$moire_data$locus)) ==
      length(missing_target_names) + length(present_target_names)
  )
  checkmate::assert(
    all(sort(present_target_names) == sort(unique(allele_freq_summary$target_name)))
  )

  one_allele_loci <- moire_object$moire_data[
    moire_object$moire_data$locus %in% missing_target_names,
  ] |>
    dplyr::select(-"sample_id") |>
    dplyr::rename(target_name = "locus", seq = "allele") |>
    dplyr::distinct() |>
    dplyr::mutate(
      post_allele_freqs_lower = 1,
      post_allele_freqs_med = 1,
      post_allele_freqs_upper = 1,
      freq = 1
    ) |>
    dplyr::select(
      "post_allele_freqs_lower",
      "post_allele_freqs_med",
      "post_allele_freqs_upper",
      "freq",
      dplyr::everything()
    )

  allele_freq_summary <- rbind(allele_freq_summary, one_allele_loci)

  target_name_count <- moire_object$moire_data |>
    dplyr::select(-"sample_id") |>
    dplyr::group_by(.data$locus) |>
    dplyr::summarise(sample_total = dplyr::n(), .groups = "drop") |>
    dplyr::rename(target_name = "locus")

  he_summary <- target_name_count |>
    dplyr::full_join(he_summary, by = "target_name")

  readr::write_tsv(coi_summary, coi_output)
  readr::write_tsv(he_summary, he_output, na = "0")
  readr::write_tsv(allele_freq_summary, allele_freq_output)
  readr::write_tsv(relatedness_summary, relatedness_output)
  readr::write_tsv(effective_coi_summary, effective_coi_output)
}

#' Assemble a named list of every estimated parameter's draws for one chain
#'
#' @param chain One element of `mcmc_results$chains`.
#' @param sample_ids Character vector of specimen IDs; order matches the
#'   per-sample draw lists (`chain$coi`, `chain$eps_pos`, etc.).
#' @param loci Character vector of locus names; order matches
#'   `chain$allele_freqs`.
#' @return A named list of numeric draw vectors, each of length
#'   `samples_per_chain`, named for the parameter it belongs to.
#' @keywords internal
extract_moire_chain_draws <- function(chain, sample_ids, loci) {
  draws <- list()
  for (s in seq_along(sample_ids)) {
    sid <- sample_ids[s]
    draws[[sprintf("coi[%s]", sid)]] <- chain$coi[[s]]
    draws[[sprintf("eps_pos[%s]", sid)]] <- chain$eps_pos[[s]]
    draws[[sprintf("eps_neg[%s]", sid)]] <- chain$eps_neg[[s]]
    # Raw (unmasked) relatedness trace, so the mixing of the sampler's
    # relatedness parameter is assessed. MOIRe masks coi <= 1 only when
    # reporting relatedness estimates, not for convergence.
    draws[[sprintf("relatedness[%s]", sid)]] <- chain$relatedness[[s]]
  }
  for (l in seq_along(loci)) {
    locus <- chain$allele_freqs[[l]]
    num_alleles <- length(locus[[1]])
    allele_freq_matrix <- matrix(unlist(locus), nrow = num_alleles)
    for (a in seq_len(num_alleles)) {
      draws[[sprintf("allele_freq[%s.%d]", loci[l], a)]] <-
        allele_freq_matrix[a, ]
    }
  }
  draws[["mean_coi"]] <- chain$mean_coi
  draws
}

#' Compute MCMC convergence diagnostics across all chains
#'
#' Assembles a posterior draws array (iteration x chain x variable) covering
#' every estimated parameter (per-sample COI, false-positive/false-negative
#' error rates, within-host relatedness; per-locus/allele frequencies; and the
#' population mean COI) and summarizes it with [summarize_convergence_draws()].
#'
#' @param mcmc_results The list returned by [run_moire()].
#' @return A data frame of convergence diagnostics, one row per parameter.
#' @keywords internal
prepare_moire_convergence_output <- function(mcmc_results) {
  check_suggested_pkg("posterior", "MOIRe convergence diagnostics")

  sample_ids <- mcmc_results$args$data$sample_ids
  loci <- mcmc_results$args$data$loci
  per_chain <- lapply(
    mcmc_results$chains,
    extract_moire_chain_draws,
    sample_ids,
    loci
  )
  var_names <- names(per_chain[[1]])
  n_iter <- length(per_chain[[1]][[1]])
  n_chains <- length(per_chain)
  draws <- array(
    NA_real_,
    dim = c(n_iter, n_chains, length(var_names)),
    dimnames = list(iteration = NULL, chain = NULL, variable = var_names)
  )
  for (i in seq_len(n_chains)) {
    draws[, i, ] <- sapply(var_names, function(v) per_chain[[i]][[v]])
  }
  summarize_convergence_draws(posterior::as_draws_array(draws))
}

#' Compute parallel tempering swap acceptance rates
#'
#' For each independent chain and each temperature rung, computes the swap
#' (exchange) acceptance rate with the adjacent hotter rung, following MOIRe's
#' own convention (see `moire::plot_chain_swaps()`):
#' `swap_acceptances / (samples_per_chain / 2)`.
#'
#' @param mcmc_results The list returned by [run_moire()].
#'
#' @details Swap acceptances are recorded per adjacent rung pair, so the rate
#'   for rung `k` describes swaps between rung `k` and rung `k + 1`; the final
#'   rung has no partner above it and its rate is `NA`. `temperature` is MOIRe's
#'   `temp_gradient` value for the rung (the power-posterior exponent in
#'   `[0, 1]`), read per chain as it may be adapted.
#'
#' @return A tibble with one row per chain-rung combination and the columns
#'   `chain`, `rung`, `temperature`, and `swap_acceptance_rate`.
#' @keywords internal
prepare_moire_acceptance_rates_output <- function(mcmc_results) {
  swap_attempts <- mcmc_results$args$samples_per_chain / 2
  chain_tables <- lapply(seq_along(mcmc_results$chains), function(chain_num) {
    chain <- mcmc_results$chains[[chain_num]]
    temps <- chain$temp_gradient
    swap_rate <- c(chain$swap_acceptances / swap_attempts, NA_real_)
    tibble::tibble(
      chain = chain_num,
      rung = seq_along(temps),
      temperature = temps,
      swap_acceptance_rate = swap_rate
    )
  })
  dplyr::bind_rows(chain_tables)
}

#' Run MOIRe from allele-table and output paths
#'
#' Reads an allele table, runs MOIRe MCMC, and writes COI, He, allele-frequency,
#' relatedness, effective-COI, and convergence summaries. Requires **moire**,
#' **checkmate**, and **posterior** (Suggests).
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele table TSV (`specimen_name`, `target_name`,
#'   `seq`). See `vignette("input-formats", package = "PGEcore")`.
#'
#' ## Outputs
#'
#' - **`coi_output`**: COI summary (`specimen_name`, `coi`, …).
#' - **`he_output`**: Heterozygosity summary (`target_name`, `he`, …).
#' - **`allele_freq_output`**: Allele frequencies (`target_name`, `seq`,
#'   `freq`, …).
#' - **`relatedness_output`**: Within-host relatedness (`specimen_name`,
#'   `within_host_rel`, …).
#' - **`effective_coi_output`**: Effective COI (`specimen_name`, `ecoi`, …).
#' - **`convergence_output`**: MCMC diagnostics (`variable`, `mean`, `median`,
#'   `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`).
#' - **`mcmc_results_output`**: Optional RDS of the full MCMC object.
#' - **`acceptance_rates_output`**: Optional PT swap rates when `pt_chains > 1`.
#'
#' ## Running
#'
#' ```r
#' moire_wrapper(
#'   allele_table = "allele_table.tsv",
#'   coi_output = "coi_output.tsv",
#'   he_output = "he_output.tsv",
#'   allele_freq_output = "allele_freq_output.tsv",
#'   relatedness_output = "relatedness_output.tsv",
#'   effective_coi_output = "effective_coi_output.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/moire_wrapper \
#'   --allele_table allele_table.tsv \
#'   --coi_output coi_output.tsv \
#'   --he_output he_output.tsv \
#'   --allele_freq_output allele_freq_output.tsv \
#'   --relatedness_output relatedness_output.tsv \
#'   --effective_coi_output effective_coi_output.tsv
#' ```
#'
#' Requires **moire**, **checkmate**, and **posterior** (Suggests).
#'
#' @param allele_table Path to allele table TSV. See *Inputs*.
#' @param allow_relatedness Logical; allow relatedness within samples.
#' @param burnin MCMC burn-in iterations.
#' @param samples_per_chain Samples per MCMC chain.
#' @param n_chains Number of independent MCMC chains. More than one is
#'   required to compute the Gelman-Rubin R-hat convergence diagnostic. This is
#'   distinct from `pt_chains` (parallel-tempering rungs within a chain).
#' @param threads Threads used to run the independent chains in
#'   parallel (MOIRe's `num_cores`).
#' @param thin Thinning interval for the MCMC sampler; only every `thin`-th
#'   sample is retained.
#' @param verbose Logical; verbose MOIRe output.
#' @param eps_pos_alpha,eps_pos_beta Positive error-rate prior.
#' @param eps_neg_alpha,eps_neg_beta Negative error-rate prior.
#' @param r_alpha,r_beta Relatedness prior.
#' @param mean_coi_shape,mean_coi_scale Mean COI prior.
#' @param max_eps_pos,max_eps_neg Maximum error rates.
#' @param record_latent_genotypes Logical; record latent genotypes.
#' @param pt_chains Number of parallel-tempering chains.
#' @param pt_grad_lower Lower bound for PT temperature gradient.
#' @param pt_num_threads Threads for parallel tempering.
#' @param adapt_temp Logical; adaptive temperature.
#' @param max_runtime Maximum MCMC runtime.
#' @param coi_output Output path for COI summary TSV. See *Outputs*.
#' @param he_output Output path for He summary TSV. See *Outputs*.
#' @param allele_freq_output Output path for allele-frequency summary TSV.
#'   See *Outputs*.
#' @param relatedness_output Output path for relatedness summary TSV. See
#'   *Outputs*.
#' @param effective_coi_output Output path for effective COI summary TSV. See
#'   *Outputs*.
#' @param mcmc_results_output Optional RDS path for full MCMC results.
#' @param convergence_output Output path for MCMC convergence diagnostics.
#'   See *Outputs*.
#' @param acceptance_rates_output Optional output path for parallel-tempering
#'   swap acceptance rates. Only meaningful when `pt_chains > 1`.
#'
#' @return Invisibly, the MOIRe MCMC result object.
#'
#' @seealso [run_moire()], `vignette("input-formats", package = "PGEcore")`
#'
#' @export
moire_wrapper <- function(allele_table,
                          allow_relatedness = TRUE,
                          burnin = 10000L,
                          samples_per_chain = 1000L,
                          n_chains = 3L,
                          threads = 1L,
                          thin = 1L,
                          verbose = FALSE,
                          eps_pos_alpha = 1,
                          eps_pos_beta = 1,
                          eps_neg_alpha = 1,
                          eps_neg_beta = 1,
                          r_alpha = 1,
                          r_beta = 1,
                          mean_coi_shape = 0.1,
                          mean_coi_scale = 10,
                          max_eps_pos = 2,
                          max_eps_neg = 2,
                          record_latent_genotypes = FALSE,
                          pt_chains = 1L,
                          pt_grad_lower = 0,
                          pt_num_threads = 1L,
                          adapt_temp = TRUE,
                          max_runtime = Inf,
                          coi_output = "coi_output.tsv",
                          he_output = "he_output.tsv",
                          allele_freq_output = "allele_freq_output.tsv",
                          relatedness_output = "relatedness_output.tsv",
                          effective_coi_output = "effective_coi_output.tsv",
                          mcmc_results_output = NULL,
                          convergence_output = "convergence_diag.tsv",
                          acceptance_rates_output = NULL) {
  check_suggested_pkg("moire", "MOIRe analysis via moire_wrapper()")
  check_suggested_pkg("checkmate", "MOIRe input validation")
  check_suggested_pkg("posterior", "MOIRe convergence diagnostics")

  if (!file.exists(allele_table)) {
    stop("allele_table file not found: ", allele_table, call. = FALSE)
  }

  moire_object <- create_moire_input(
    allele_table,
    allow_relatedness,
    burnin,
    samples_per_chain,
    thin,
    verbose,
    eps_pos_alpha,
    eps_pos_beta,
    eps_neg_alpha,
    eps_neg_beta,
    r_alpha,
    r_beta,
    mean_coi_shape,
    mean_coi_scale,
    max_eps_pos,
    max_eps_neg,
    record_latent_genotypes,
    pt_chains,
    pt_grad_lower,
    pt_num_threads,
    adapt_temp,
    max_runtime,
    n_chains,
    threads
  )

  moire_results <- run_moire(moire_object)
  summarize_and_write_moire_results(
    moire_object,
    moire_results,
    coi_output,
    he_output,
    allele_freq_output,
    relatedness_output,
    effective_coi_output
  )

  readr::write_tsv(
    prepare_moire_convergence_output(moire_results),
    convergence_output
  )

  if (!is.null(acceptance_rates_output)) {
    if (length(moire_results$chains[[1]]$temp_gradient) > 1) {
      readr::write_tsv(
        prepare_moire_acceptance_rates_output(moire_results),
        acceptance_rates_output
      )
    } else {
      warning(
        "acceptance_rates_output was provided but parallel tempering was not ",
        "used (pt_chains must be > 1); no acceptance rates table written.",
        call. = FALSE
      )
    }
  }

  if (!is.null(mcmc_results_output)) {
    saveRDS(moire_results, mcmc_results_output)
  }

  invisible(moire_results)
}
