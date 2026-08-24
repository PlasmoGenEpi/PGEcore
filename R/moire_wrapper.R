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
                               max_runtime) {
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
  if (pt_chains > 1) {
    pt_chains <- seq(from = pt_grad_lower, to = 1, length.out = pt_chains)
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

  p <- moire_object$moire_parameters
  checkmate::assert_logical(p$allow_relatedness, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$burnin, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$samples_per_chain, any.missing = FALSE, len = 1)
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
  checkmate::assert_numeric(p$pt_num_threads, any.missing = FALSE, len = 1)
  checkmate::assert_logical(p$adapt_temp, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(p$max_runtime, any.missing = FALSE, len = 1)

  message("Returning Moire object")
  moire_object
}

#' Run MOIRe MCMC analysis
#'
#' The **moire** package is an optional dependency (Suggests). It is not
#' installed automatically with PGEcore.
#'
#' @param moire_object List created by [create_moire_input()] (or
#'   [moire_wrapper()]), with `moire_data` and `moire_parameters`.
#'
#' @return The object returned by [moire::run_mcmc()].
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
}

#' Summarize MOIRe MCMC results and write TSV files
#'
#' @keywords internal
summarize_and_write_moire_results <- function(moire_object,
                                              mcmc_results,
                                              coi_summary_o,
                                              he_summary_o,
                                              allele_freq_summary_o,
                                              relatedness_summary_o,
                                              effective_coi_summary_o) {
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

  readr::write_tsv(coi_summary, coi_summary_o)
  readr::write_tsv(he_summary, he_summary_o, na = "0")
  readr::write_tsv(allele_freq_summary, allele_freq_summary_o)
  readr::write_tsv(relatedness_summary, relatedness_summary_o)
  readr::write_tsv(effective_coi_summary, effective_coi_summary_o)
}

#' Run MOIRe from allele-table and output paths
#'
#' File-oriented entry point used by the `moire_wrapper` CLI. Optional
#' **moire** (and **checkmate** for input checks) must be installed separately.
#'
#' @param allele_table Path to a TSV with columns `specimen_name`,
#'   `target_name`, and `seq`.
#' @param allow_relatedness Logical; allow relatedness within samples.
#' @param burnin MCMC burn-in iterations.
#' @param samples_per_chain Samples per MCMC chain.
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
#' @param coi_summary Output path for COI summary TSV.
#' @param he_summary Output path for He summary TSV.
#' @param allele_freq_summary Output path for allele-frequency summary TSV.
#' @param relatedness_summary Output path for relatedness summary TSV.
#' @param effective_coi_summary Output path for effective COI summary TSV.
#' @param mcmc_results_output Optional RDS path for full MCMC results.
#'
#' @return Invisibly, the MOIRe MCMC result object.
#' @export
moire_wrapper <- function(allele_table,
                          allow_relatedness = TRUE,
                          burnin = 10000L,
                          samples_per_chain = 1000L,
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
                          coi_summary = "coi_summary.tsv",
                          he_summary = "he_summary.tsv",
                          allele_freq_summary = "allele_freq_summary.tsv",
                          relatedness_summary = "relatedness_summary.tsv",
                          effective_coi_summary = "effective_coi_summary.tsv",
                          mcmc_results_output = NULL) {
  check_suggested_pkg("moire", "MOIRe analysis via moire_wrapper()")
  check_suggested_pkg("checkmate", "MOIRe input validation")

  if (!file.exists(allele_table)) {
    stop("allele_table file not found: ", allele_table, call. = FALSE)
  }

  moire_object <- create_moire_input(
    allele_table,
    allow_relatedness,
    burnin,
    samples_per_chain,
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
    max_runtime
  )

  moire_results <- run_moire(moire_object)
  summarize_and_write_moire_results(
    moire_object,
    moire_results,
    coi_summary,
    he_summary,
    allele_freq_summary,
    relatedness_summary,
    effective_coi_summary
  )

  if (!is.null(mcmc_results_output)) {
    saveRDS(moire_results, mcmc_results_output)
  }

  invisible(moire_results)
}
