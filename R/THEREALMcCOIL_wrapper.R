#' Read SNP call data for THEREALMcCOIL
#'
#' @param snp_calls Path to an independent, collapsed SNP call TSV.
#' @return A data frame with `specimen_name`, `snp_name`, `seq_base`, `reads`.
#' @keywords internal
read_and_preprocess_snp_call <- function(snp_calls) {
  required_cols <- c("specimen_name", "snp_name", "reads", "seq_base")
  df_snp_call <- readr::read_tsv(
    snp_calls,
    col_types = readr::cols(specimen_name = readr::col_character())
  )
  validate_required_columns(df_snp_call, required_cols, "SNP data")
  df_snp_call |>
    dplyr::select(dplyr::all_of(required_cols))
}

#' Format SNP calls for the McCOIL categorical model
#'
#' @param df Output of `read_and_preprocess_snp_call()`.
#' @return A data frame (samples x sites) of scores `1` / `0` / `0.5` / `-1`.
#' @keywords internal
prep_input_categorical <- function(df) {
  major_allele <- df |>
    dplyr::group_by(.data$snp_name, .data$seq_base) |>
    dplyr::summarize(total = sum(.data$reads), .groups = "keep") |>
    dplyr::group_by(.data$snp_name) |>
    dplyr::slice_max(.data$total, n = 1, with_ties = FALSE) |>
    dplyr::select("snp_name", "seq_base") |>
    dplyr::rename(major_allele = "seq_base")

  df_recode <- df |>
    dplyr::left_join(major_allele, by = "snp_name") |>
    dplyr::mutate(allele_idx = dplyr::if_else(.data$seq_base == .data$major_allele, 1, 0)) |>
    dplyr::select("specimen_name", "snp_name", "allele_idx") |>
    dplyr::group_by(.data$specimen_name, .data$snp_name) |>
    dplyr::summarise(
      count = dplyr::n(),
      score = dplyr::case_when(
        .data$count == 0 ~ -1,
        all(.data$allele_idx == 0) ~ 0,
        all(.data$allele_idx == 1) ~ 1,
        TRUE ~ 0.5
      ),
      .groups = "drop"
    )

  df_wide <- df_recode |>
    dplyr::select(-"count") |>
    tidyr::pivot_wider(names_from = "snp_name", values_from = "score")

  df_wide[is.na(df_wide)] <- -1

  df_mat <- data.frame(df_wide[, -1])
  colnames(df_mat) <- colnames(df_wide)[-1]
  rownames(df_mat) <- dplyr::pull(df_wide, "specimen_name")
  df_mat
}

#' Restrict SNP calls to biallelic loci for the proportional model
#'
#' The proportional model represents each locus with exactly two alleles
#' (`a1`/`a2`), so loci with any other number of alleles cannot be encoded
#' correctly. Loci with exactly two distinct `seq_base` values are kept;
#' monomorphic and multiallelic loci are dropped with a warning.
#'
#' @param df Output of `read_and_preprocess_snp_call()`.
#' @return `df` filtered to biallelic loci only.
#' @keywords internal
filter_biallelic <- function(df) {
  allele_counts <- df |>
    dplyr::distinct(.data$snp_name, .data$seq_base) |>
    dplyr::count(.data$snp_name, name = "n_alleles")

  n_mono <- sum(allele_counts$n_alleles < 2)
  n_multi <- sum(allele_counts$n_alleles > 2)
  if (n_mono > 0 || n_multi > 0) {
    warning(
      "Dropping non-biallelic loci for the proportional model: ",
      n_mono, " monomorphic, ", n_multi, " multiallelic. ",
      "The proportional model supports only biallelic loci.",
      call. = FALSE
    )
  }

  biallelic_loci <- allele_counts |>
    dplyr::filter(.data$n_alleles == 2) |>
    dplyr::pull("snp_name")
  if (length(biallelic_loci) == 0) {
    stop(
      "No biallelic loci remain after filtering; the proportional model ",
      "cannot be run on this input.",
      call. = FALSE
    )
  }

  dplyr::filter(df, .data$snp_name %in% biallelic_loci)
}

#' Format SNP calls for the McCOIL proportional model
#'
#' @param df Output of `read_and_preprocess_snp_call()`.
#' @return A list with `a1` and `a2` read-count matrices.
#' @keywords internal
prep_input_prop <- function(df) {
  df <- filter_biallelic(df)

  # Assign the two alleles of each locus to index 1 or 2 within that locus, so
  # the assignment cannot depend on the number of alleles seen at other loci.
  allele_map <- df |>
    dplyr::distinct(.data$snp_name, .data$seq_base) |>
    dplyr::arrange(.data$snp_name, .data$seq_base) |>
    dplyr::group_by(.data$snp_name) |>
    dplyr::mutate(allele_idx = dplyr::row_number()) |>
    dplyr::ungroup()
  df_with_allele_idx <- df |>
    dplyr::left_join(allele_map, by = c("snp_name", "seq_base"))

  df_allele1 <- df_with_allele_idx |>
    dplyr::filter(.data$allele_idx == 1) |>
    dplyr::select("specimen_name", "snp_name", "reads") |>
    tidyr::pivot_wider(values_from = "reads", names_from = "snp_name") |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., 0)))

  column_names <- colnames(df_allele1)
  df_allele1 <- data.frame(df_allele1)
  colnames(df_allele1) <- column_names
  row.names(df_allele1) <- df_allele1$specimen_name
  df_allele1 <- df_allele1[, -1, drop = FALSE]

  df_allele2 <- df_with_allele_idx |>
    dplyr::filter(.data$allele_idx == 2) |>
    dplyr::select("specimen_name", "snp_name", "reads") |>
    tidyr::pivot_wider(values_from = "reads", names_from = "snp_name") |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., 0)))

  column_names <- colnames(df_allele2)
  df_allele2 <- data.frame(df_allele2)
  colnames(df_allele2) <- column_names
  row.names(df_allele2) <- df_allele2$specimen_name
  df_allele2 <- df_allele2[, -1, drop = FALSE]

  # The two subsets are pivoted independently, so each names its columns in the
  # order it happens to meet the loci and carries only the specimens that have a
  # read for that allele. McCOIL_prop sizes both matrices from a1 and reads them
  # as flat vectors, so anything not on a shared grid pairs one locus's allele-1
  # count with another locus's allele-2 count -- the likelihood never improves
  # and every specimen stays at its starting COI. Put both on the full
  # specimen x locus grid; an absent combination is a count of zero, and no
  # specimen is dropped for lacking one of the two alleles.
  specimens <- unique(df_with_allele_idx$specimen_name)
  loci <- unique(df_with_allele_idx$snp_name)
  on_grid <- function(d) {
    out <- matrix(0, nrow = length(specimens), ncol = length(loci),
                  dimnames = list(specimens, loci))
    r <- intersect(specimens, rownames(d))
    cl <- intersect(loci, colnames(d))
    if (length(r) > 0 && length(cl) > 0) {
      out[r, cl] <- as.matrix(d[r, cl, drop = FALSE])
    }
    as.data.frame(out)
  }

  a1 <- on_grid(df_allele1)
  a2 <- on_grid(df_allele2)

  # McCOIL_prop reads a locus as missing only when a count is negative; its
  # likelihood divides by the two counts, so a locus with no reads for either
  # allele gives 0/0 and turns the whole per-specimen sum into NaN. No proposal
  # then clears the acceptance test and that specimen keeps its starting COI for
  # the entire chain. Mark uncovered specimen-locus pairs as missing instead.
  uncovered <- (a1 + a2) == 0
  a1[uncovered] <- -1
  a2[uncovered] <- -1

  list(a1 = a1, a2 = a2)
}

#' Run a single McCOIL chain and return its per-iteration trace
#'
#' Runs one MCMC chain of the requested model with a given seed, writing its
#' trace and summary to `output` under `work_dir`, and reads the trace back.
#'
#' @param mccoil_input Prepared model input: the genotype matrix for the
#'   categorical model, or a list with `a1`/`a2` matrices for the proportional
#'   model.
#' @param model `"categorical"` or `"proportional"`.
#' @param seed Random seed for this chain.
#' @param work_dir Directory for McCOIL temp traces.
#' @param output Trace filename (under `work_dir`) for this chain.
#' @inheritParams THEREALMcCOIL_wrapper
#'
#' @return A data frame of the raw per-iteration trace, one row per iteration
#'   plus the trailing acceptance-count row.
#' @keywords internal
run_mccoil_chain <- function(mccoil_input,
                             model,
                             maxCOI,
                             threshold_ind,
                             threshold_site,
                             totalrun,
                             burnin,
                             M0,
                             e1,
                             e2,
                             epsilon,
                             err_method,
                             seed,
                             work_dir,
                             output) {
  # The compiled McCOIL code draws from R's RNG stream via GetRNGstate(), so
  # seeding immediately before the .C() call fixes the chain reproducibly.
  set.seed(seed)
  if (model == "categorical") {
    run_mccoil_categorical(
      mccoil_input,
      maxCOI = maxCOI,
      threshold_ind = threshold_ind,
      threshold_site = threshold_site,
      totalrun = totalrun,
      burnin = burnin,
      M0 = M0,
      e1 = e1,
      e2 = e2,
      err_method = err_method,
      path = work_dir,
      output = output
    )
  } else {
    run_mccoil_proportional(
      mccoil_input$a1,
      mccoil_input$a2,
      maxCOI = maxCOI,
      totalrun = totalrun,
      burnin = burnin,
      M0 = M0,
      epsilon = epsilon,
      err_method = err_method,
      path = work_dir,
      output = output
    )
  }
  utils::read.table(file.path(work_dir, output), header = FALSE)
}

#' Run several independent McCOIL chains into `work_dir`
#'
#' Prepares the model input once, then runs `n_chains` chains with seeds
#' `seed, seed + 1, ...`. Chains are independent and write to distinct trace
#' files in `work_dir`, so they can run in parallel forks; on Windows, where
#' forking is unavailable, they run sequentially.
#'
#' @param df Preprocessed SNP calls.
#' @param model `"categorical"` or `"proportional"`.
#' @param work_dir Directory for McCOIL temp traces (typically under
#'   `tempdir()`).
#' @param output Base filename for chain 1, written under `work_dir`; later
#'   chains append a `_chain<i>` suffix.
#' @inheritParams THEREALMcCOIL_wrapper
#'
#' @return A list with `traces` (per-chain trace data frames, chain 1 first)
#'   and `summary_path`, the summary file written by chain 1.
#' @keywords internal
run_mccoil_chains <- function(df,
                              model = "categorical",
                              maxCOI = 25,
                              threshold_ind = 20,
                              threshold_site = 20,
                              totalrun = 10000,
                              burnin = 1000,
                              M0 = 15,
                              e1 = 0.05,
                              e2 = 0.05,
                              epsilon = 0.02,
                              err_method = 1,
                              seed = 321,
                              n_chains = 3,
                              work_dir,
                              output = "McCOIL_out.txt") {
  if (!model %in% c("categorical", "proportional")) {
    stop("--model must be one of categorical|proportional", call. = FALSE)
  }
  if (!err_method %in% c(1, 3)) {
    stop("--err_method must be one of 1|3", call. = FALSE)
  }
  n_chains <- as.integer(n_chains)
  if (is.na(n_chains) || n_chains < 1L) {
    stop("--n_chains must be a positive integer", call. = FALSE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  mccoil_input <- if (model == "categorical") {
    prep_input_categorical(df)
  } else {
    prep_input_prop(df)
  }

  # Chain 1 keeps the base filename, so its summary drives the COI/SLAF output.
  chains <- seq_len(n_chains)
  output_names <- ifelse(
    chains == 1L,
    output,
    paste0(sub("\\.txt$", "", output), "_chain", chains, ".txt")
  )
  seeds <- seed + chains - 1

  run_chain <- function(i) {
    run_mccoil_chain(
      mccoil_input,
      model = model,
      maxCOI = maxCOI,
      threshold_ind = threshold_ind,
      threshold_site = threshold_site,
      totalrun = totalrun,
      burnin = burnin,
      M0 = M0,
      e1 = e1,
      e2 = e2,
      epsilon = epsilon,
      err_method = err_method,
      seed = seeds[i],
      work_dir = work_dir,
      output = output_names[i]
    )
  }

  traces <- if (n_chains > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(chains, run_chain, mc.cores = n_chains)
  } else {
    lapply(chains, run_chain)
  }

  failed <- which(!vapply(traces, is.data.frame, logical(1)))
  if (length(failed) > 0) {
    condition <- attr(traces[[failed[1]]], "condition")
    stop(
      "McCOIL chain(s) ", paste(failed, collapse = ", "), " failed: ",
      if (is.null(condition)) "unknown error" else conditionMessage(condition),
      call. = FALSE
    )
  }

  list(
    traces = traces,
    summary_path = file.path(work_dir, paste0(output, "_summary.txt"))
  )
}

#' Compute MCMC convergence diagnostics across McCOIL chains
#'
#' Assembles a posterior draws array from the post-burn-in portion of the
#' per-chain traces and summarises it with [summarize_convergence_draws()].
#'
#' @param traces Per-chain trace data frames from [run_mccoil_chains()].
#' @param summary_path Path to the `*_summary.txt` written by chain 1, used for
#'   the parameter names and their order in the traces.
#' @inheritParams THEREALMcCOIL_wrapper
#'
#' @return A data frame of convergence diagnostics, one row per parameter.
#' @keywords internal
prepare_mccoil_convergence_output <- function(traces, summary_path, totalrun, burnin) {
  check_suggested_pkg("posterior", "MCMC convergence diagnostics")
  summary_df <- utils::read.table(
    summary_path,
    sep = "\t",
    header = TRUE,
    colClasses = c(name = "character")
  )
  coi_names <- summary_df$name[summary_df$CorP == "C"]
  freq_names <- summary_df$name[summary_df$CorP == "P"]
  err_names <- summary_df$name[!summary_df$CorP %in% c("C", "P")]
  n <- length(coi_names)
  k <- length(freq_names)
  n_err <- length(err_names)

  # Trace layout: column 1 is the iteration index; columns 2:(n + 1) are COI
  # per specimen; the next k are allele frequency per locus; any remaining
  # columns are error parameters. Rows 1:totalrun are iterations (the trailing
  # acceptance-count row is excluded); the first `burnin` are dropped.
  keep_rows <- (burnin + 1):totalrun
  keep_cols <- 2:(1 + n + k + n_err)
  var_names <- c(
    paste0("coi[", coi_names, "]"),
    paste0("freq[", freq_names, "]"),
    err_names
  )

  draws <- array(
    NA_real_,
    dim = c(length(keep_rows), length(traces), length(var_names)),
    dimnames = list(iteration = NULL, chain = NULL, variable = var_names)
  )
  for (i in seq_along(traces)) {
    draws[, i, ] <- as.matrix(traces[[i]][keep_rows, keep_cols])
  }
  summarize_convergence_draws(posterior::as_draws_array(draws))
}

#' Format McCOIL summary TSV into PGE COI and SLAF tables
#'
#' @param summary_path Path to `*_summary.txt` written by McCOIL.
#' @return A list with `slaf` and `coi` tibbles.
#' @keywords internal
format_mccoil_output <- function(summary_path) {
  df_mccoil <- utils::read.table(
    summary_path,
    sep = "\t",
    header = TRUE,
    colClasses = c(name = "character")
  )

  df_slaf <- df_mccoil |>
    dplyr::filter(.data$CorP == "P") |>
    dplyr::select("name", "median") |>
    dplyr::rename(variant = "name", freq = "median")

  df_coi <- df_mccoil |>
    dplyr::filter(.data$CorP == "C") |>
    dplyr::select("name", "median") |>
    dplyr::rename(specimen_name = "name", coi = "median")

  list(slaf = df_slaf, coi = df_coi)
}

#' Write formatted McCOIL output
#'
#' @param df_formated List from `format_mccoil_output()`.
#' @param slaf_path SLAF TSV path (`variant`, `freq`).
#' @param coi_path COI TSV path (`specimen_name`, `coi`).
#' @keywords internal
write_mccoil_output <- function(df_formated, slaf_path, coi_path) {
  readr::write_tsv(df_formated$slaf, slaf_path)
  readr::write_tsv(df_formated$coi, coi_path)
}

#' Remove McCOIL intermediate files from a working directory
#'
#' @param work_dir Directory that may contain `McCOIL*` traces.
#' @keywords internal
clean_up_mccoil <- function(work_dir) {
  if (!dir.exists(work_dir)) {
    return(invisible(NULL))
  }
  traces <- list.files(work_dir, pattern = "^McCOIL", full.names = TRUE)
  if (length(traces) > 0) {
    invisible(file.remove(traces))
  }
  invisible(NULL)
}

mccoil_blank <- function(x) {
  missing(x) || is.null(x) || (is.character(x) && !nzchar(x))
}

#' Estimate COI and allele frequencies with THEREALMcCOIL
#'
#' Runs THEREALMcCOIL MCMC on SNP calls to estimate per-specimen COI and
#' single-locus allele frequencies. Compiled C routines (`McCOIL_categorical`,
#' `McCOIL_prop`) are linked at package install time. Requires **posterior**
#' (Suggests) for convergence diagnostics.
#'
#' ## Inputs
#'
#' - **`snp_calls`**: SNP-calls TSV (at least `specimen_name`, `snp_name`,
#'   `reads`, `seq_base`; typically also `target_name`, `pos`, `he`). See
#'   `vignette("input-formats", package = "PGEcore")`.
#'
#' ## Outputs
#'
#' - **`slaf_output`**: Single-locus allele frequencies (`variant`, `freq`).
#' - **`coi_output`**: COI estimates (`specimen_name`, `coi`).
#' - **`convergence_output`**: MCMC diagnostics (`variable`, `mean`, `median`,
#'   `sd`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`).
#'
#' ## Running
#'
#' ```r
#' THEREALMcCOIL_wrapper(
#'   snp_calls = "snp_calls.tsv",
#'   slaf_output = "slaf.tsv",
#'   coi_output = "coi.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/THEREALMcCOIL_wrapper \
#'   --snp_calls snp_calls.tsv \
#'   --slaf_output slaf.tsv \
#'   --coi_output coi.tsv
#' ```
#'
#' @param snp_calls Path to SNP-calls TSV. See *Inputs*.
#' @param slaf_output Output TSV of allele frequencies. See *Outputs*.
#' @param coi_output Output TSV of COI estimates. See *Outputs*.
#' @param model `"categorical"` (heterozygous/homozygous calls) or
#'   `"proportional"` (allele frequency / read-count data).
#' @param maxCOI Upper bound for COI.
#' @param threshold_ind Minimum sites per sample (categorical model).
#' @param threshold_site Minimum samples per locus (categorical model).
#' @param totalrun Total MCMC iterations.
#' @param burnin Burn-in iterations.
#' @param M0 Initial COI.
#' @param e1 Probability of calling homozygous loci heterozygous (categorical).
#' @param e2 Probability of calling heterozygous loci homozygous (categorical).
#' @param epsilon Error parameter for the proportional model.
#' @param err_method `1`: treat error rates as constants; `3`: estimate them
#'   with COI and allele frequencies.
#' @param seed Random seed for the first chain; chain `i` uses `seed + i - 1`.
#' @param n_chains Number of independent MCMC chains. More than one chain is
#'   required for the Gelman-Rubin R-hat diagnostic.
#' @param convergence_summary_output Output TSV of the run-level convergence
#'   summary. See *Outputs*.
#' @param convergence_output Output TSV of MCMC convergence diagnostics. See
#'   *Outputs*.
#'
#' @return A list with `slaf`, `coi` and `convergence` (invisibly after writing
#'   outputs).
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
THEREALMcCOIL_wrapper <- function(snp_calls,
                                  slaf_output,
                                  coi_output,
                                  model = "categorical",
                                  maxCOI = 25L,
                                  threshold_ind = 20L,
                                  threshold_site = 20L,
                                  totalrun = 10000L,
                                  burnin = 1000L,
                                  M0 = 15L,
                                  e1 = 0.05,
                                  e2 = 0.05,
                                  epsilon = 0.02,
                                  err_method = 1L,
                                  seed = 321L,
                                  n_chains = 3L,
                                  convergence_output = "convergence_diag.tsv",
                                  convergence_summary_output =
                                    "convergence_summary.tsv") {
  if (mccoil_blank(snp_calls)) {
    stop("--snp_calls must be set", call. = FALSE)
  }
  if (mccoil_blank(slaf_output)) {
    stop("--slaf_output must be set", call. = FALSE)
  }
  if (mccoil_blank(coi_output)) {
    stop("--coi_output must be set", call. = FALSE)
  }
  if (mccoil_blank(convergence_output)) {
    stop("--convergence_output must be set", call. = FALSE)
  }
  # Fail before the MCMC runs rather than after, when diagnostics are written.
  check_suggested_pkg("posterior", "MCMC convergence diagnostics")

  df <- read_and_preprocess_snp_call(snp_calls)
  work_dir <- tempfile("McCOIL_")
  dir.create(work_dir)
  on.exit(clean_up_mccoil(work_dir), add = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  chains <- run_mccoil_chains(
    df,
    model = model,
    maxCOI = maxCOI,
    threshold_ind = threshold_ind,
    threshold_site = threshold_site,
    totalrun = totalrun,
    burnin = burnin,
    M0 = M0,
    e1 = e1,
    e2 = e2,
    epsilon = epsilon,
    err_method = err_method,
    seed = seed,
    n_chains = n_chains,
    work_dir = work_dir,
    output = "McCOIL_out.txt"
  )
  df_formated <- format_mccoil_output(chains$summary_path)
  write_mccoil_output(df_formated, slaf_output, coi_output)

  convergence <- prepare_mccoil_convergence_output(
    chains$traces,
    chains$summary_path,
    totalrun = totalrun,
    burnin = burnin
  )
  readr::write_tsv(convergence, convergence_output)
  convergence_summary <- summarize_convergence_run(convergence)
  readr::write_tsv(convergence_summary, convergence_summary_output)

  invisible(c(df_formated, list(
    convergence = convergence, convergence_summary = convergence_summary
  )))
}
