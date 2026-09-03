#' Read an allele table with counts for SNP-Slice
#'
#' @keywords internal
create_snpslice_allele_table_input <- function(allele_table_path,
                                               specimen_name_col = "specimen_name",
                                               target_name_col = "target_name",
                                               target_value_col = "aa",
                                               target_count_col = "reads") {
  allele_table <- readr::read_tsv(
    allele_table_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      !!target_count_col := readr::col_double()
    ),
    progress = FALSE
  ) |>
    dplyr::select(dplyr::all_of(c(
      specimen_name_col,
      target_name_col,
      target_value_col,
      target_count_col
    ))) |>
    dplyr::rename(
      specimen_name = dplyr::all_of(specimen_name_col),
      target_name = dplyr::all_of(target_name_col),
      target_value = dplyr::all_of(target_value_col),
      target_count = dplyr::all_of(target_count_col)
    )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(target_value),
    is.double(target_count),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(target_value),
    !is.na(target_count)
  )
  stop_on_validate_fails(allele_table, rules, "allele_table")
  allele_table
}

#' Read loci groups for SNP-Slice multilocus frequencies
#'
#' @keywords internal
create_snpslice_loci_group_input <- function(loci_groups_path,
                                             allele_table,
                                             target_name_col = "target_name") {
  stopifnot(is.character(loci_groups_path))

  loci_groups <- readr::read_tsv(
    loci_groups_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      aa_position = readr::col_integer()
    ),
    progress = FALSE
  ) |>
    dplyr::select("group_id", dplyr::all_of(target_name_col)) |>
    dplyr::rename(target_name = dplyr::all_of(target_name_col))

  rules <- validate::validator(
    is.character(group_id),
    is.character(target_name),
    !is.na(group_id),
    !is.na(target_name)
  )
  stop_on_validate_fails(loci_groups, rules, "loci_groups")

  loci_groups <- split(loci_groups$target_name, loci_groups$group_id)

  for (lg in names(loci_groups)) {
    missing_trgs <- setdiff(loci_groups[[lg]], allele_table$target_name)
    if (length(missing_trgs) > 0) {
      warning(
        "The target(s) ",
        paste(missing_trgs, collapse = ", "),
        " in the group ",
        lg,
        " are missing and will be excluded.",
        call. = FALSE
      )
    }
    non_biallelic_trgs <- allele_table |>
      dplyr::filter(.data$target_name %in% loci_groups[[lg]]) |>
      dplyr::group_by(.data$target_name) |>
      dplyr::filter(dplyr::n_distinct(.data$target_value) > 2) |>
      dplyr::pull("target_name") |>
      unique()
    if (length(non_biallelic_trgs) > 0) {
      warning(
        "The target(s) ",
        paste(non_biallelic_trgs, collapse = ", "),
        " in the group ",
        lg,
        " have more than two alleles and will be excluded.",
        call. = FALSE
      )
    }
    loci_groups[[lg]] <- setdiff(
      loci_groups[[lg]],
      c(missing_trgs, non_biallelic_trgs)
    )
  }
  if (length(unlist(loci_groups)) == 0) {
    stop("The data is missing all loci in loci groups", call. = FALSE)
  }
  loci_groups
}

#' Format SNP-Slice allele frequencies with variantstring names
#'
#' @keywords internal
prepare_snpslice_af_output <- function(snpslice_res, loci_groups, estimator) {
  format_af_table_w_variantstring <- function(af_table, group_id, loci_groups) {
    prep_variantstring_input <- function(allele, loci_names) {
      aa <- stringr::str_split_1(allele, "\\|")
      tibble::tibble(gene_pos = loci_names, aa = aa) |>
        tidyr::separate_wider_delim(
          "gene_pos",
          ":",
          names = c("gene", "pos")
        ) |>
        dplyr::mutate(pos = as.integer(.data$pos)) |>
        dplyr::mutate(
          n_aa = 1,
          het = FALSE,
          phased = TRUE,
          read_count = NA
        ) |>
        dplyr::relocate("aa", .before = "read_count")
    }

    tibble::as_tibble(af_table) |>
      dplyr::filter(.data$frequency > 0) |>
      dplyr::mutate(
        allele = lapply(
          .data$allele,
          prep_variantstring_input,
          loci_groups[[group_id]]
        )
      ) |>
      dplyr::mutate(allele = variantstring::long_to_variant(.data$allele))
  }

  snp_slicer_af_by_group <- snp.slicer::calculate_allele_frequencies_by_sets(
    snpslice_res,
    loci_groups,
    estimate = estimator
  )
  af_tib <- tibble::tibble(
    group_id = names(snp_slicer_af_by_group),
    af_tib = unname(snp_slicer_af_by_group)
  ) |>
    dplyr::mutate(
      af_tib = purrr::map2(
        .data$af_tib,
        .data$group_id,
        format_af_table_w_variantstring,
        loci_groups
      )
    ) |>
    tidyr::unnest("af_tib") |>
    dplyr::rename(variant = "allele", freq = "frequency")

  if (identical(estimator, "posterior")) {
    af_tib <- dplyr::select(af_tib, -"mean_count", -"n_samples")
  } else {
    af_tib <- dplyr::rename(
      af_tib,
      allele_total = "total_parasites",
      allele_count = "count"
    )
  }
  af_tib
}

#' Consensus COI across SNP-Slice restarts
#'
#' Restarts are independent optimisations that land on different strain
#' dictionaries, so a strain index means nothing across them. Haplotypes are
#' matched by their dictionary row instead, and duplicate rows within a restart
#' are collapsed before counting. `membership[i, h]` is then the fraction of
#' restarts in which host `i` carries haplotype `h`.
#'
#' Consensus support is `colSums(membership)`, and each haplotype is weighted by
#' `1 - exp(-support / mean support)`. Support scales with how finely the loci
#' resolve strains, so the weight is normalised by the mean rather than by a
#' fixed count.
#'
#' @param chains List of per-restart results, or a single result object.
#' @param estimate Point estimate to read from each restart, `"map"` or
#'   `"final_sample"`.
#' @return Numeric vector, one consensus COI per host, floored at 1.
#' @keywords internal
snpslice_consensus_coi <- function(chains, estimate) {
  if (!is.null(chains) && !is.null(chains$map_allocation_matrix)) {
    chains <- list(chains)
  }
  per_chain <- lapply(chains, function(ch) {
    mats <- snp.slicer:::point_estimate_matrices(ch, estimate)
    A <- mats$A
    D <- mats$D
    haplotype <- apply(D, 1L, paste0, collapse = "")
    columns <- split(seq_len(ncol(A)), haplotype)
    lapply(columns, function(k) which(rowSums(A[, k, drop = FALSE]) > 0))
  })
  n_hosts <- nrow(snp.slicer:::point_estimate_matrices(chains[[1L]], estimate)$A)
  haplotypes <- unique(unlist(lapply(per_chain, names), use.names = FALSE))
  membership <- matrix(0, nrow = n_hosts, ncol = length(haplotypes),
                       dimnames = list(NULL, haplotypes))
  for (chain in per_chain) {
    for (h in names(chain)) {
      membership[chain[[h]], h] <- membership[chain[[h]], h] + 1
    }
  }
  membership <- membership / length(per_chain)
  support <- colSums(membership)
  mean_support <- if (length(support) == 0L) 0 else mean(support)
  weights <- if (mean_support > 0) {
    1 - exp(-support / mean_support)
  } else {
    rep(0, length(support))
  }
  pmax(as.vector(membership %*% weights), 1)
}

#' Format SNP-Slice COI estimates
#'
#' @keywords internal
prepare_snpslice_coi_output <- function(snpslice_res,
                                        specimen_name_col,
                                        estimator) {
  coi_tib <- snp.slicer::calculate_individual_coi(
    snpslice_res,
    estimate = estimator
  ) |>
    dplyr::select(-"host_index") |>
    dplyr::rename(!!specimen_name_col := "host_id", coi = "coi_estimate")
  if (!identical(estimator, "posterior")) {
    coi_tib <- dplyr::select(coi_tib, -"coi_sd", -"coi_lower", -"coi_upper")
  }
  chains <- if (is.null(snpslice_res$chains)) {
    snp.slicer::get_chain(snpslice_res, NULL)
  } else {
    snpslice_res$chains
  }
  consensus <- snpslice_consensus_coi(
    chains,
    if (identical(estimator, "posterior")) "map" else estimator
  )
  if (length(consensus) != nrow(coi_tib)) {
    stop("Consensus COI length does not match the per-host COI table.",
         call. = FALSE)
  }
  coi_tib$coi_cons_weighted <- consensus
  dplyr::relocate(coi_tib, "coi_cons_weighted", .after = "coi")
}

#' Lin's concordance correlation coefficient
#'
#' Returns `NA_real_` for fewer than three complete pairs and 1 when both
#' vectors are constant and equal.
#'
#' @keywords internal
snpslice_ccc <- function(x, y) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 3L) {
    return(NA_real_)
  }
  vx <- stats::var(x) * (n - 1) / n
  vy <- stats::var(y) * (n - 1) / n
  cxy <- stats::cov(x, y) * (n - 1) / n
  denom <- vx + vy + (mean(x) - mean(y))^2
  if (denom == 0) {
    return(1)
  }
  2 * cxy / denom
}

#' Format SNP-Slice per-restart optimisation diagnostics
#'
#' One row per restart: `chain_id`, `seed`, `map_logpost`, `is_best`,
#' `map_iteration`, `final_iteration`, `plateau_frac` (`map_iteration` divided
#' by `final_iteration`), `map_kstar`, `map_ktrunc`, `coi_mean`, and
#' `coi_ccc_to_best` (Lin's CCC between that restart's per-host COI and the
#' reported restart's).
#'
#' @keywords internal
prepare_snpslice_optim_output <- function(snpslice_res) {
  # snp_slice() returns the reported chain at the top level and only attaches
  # $chains when n_chains > 1.
  chains <- snpslice_res$chains
  if (is.null(chains)) {
    chains <- list(snpslice_res)
  }
  best <- snpslice_res$best_chain
  if (is.null(best)) {
    best <- 1L
  }
  # COI per host is the row sum of the MAP allocation matrix.
  coi_by_chain <- lapply(chains, function(ch) rowSums(ch$map_allocation_matrix))
  coi_best <- coi_by_chain[[best]]

  purrr::list_rbind(purrr::imap(chains, function(ch, i) {
    d <- ch$diagnostics
    tibble::tibble(
      chain_id = if (is.null(d$chain_id)) i else d$chain_id,
      seed = if (is.null(d$seed)) NA_integer_ else d$seed,
      map_logpost = d$map_logpost,
      is_best = identical(as.integer(i), as.integer(best)),
      map_iteration = d$map_iteration,
      final_iteration = d$final_iteration,
      plateau_frac = d$map_iteration / d$final_iteration,
      map_kstar = d$map_kstar,
      map_ktrunc = d$map_ktrunc,
      coi_mean = mean(coi_by_chain[[i]]),
      coi_ccc_to_best = snpslice_ccc(coi_by_chain[[i]], coi_best)
    )
  }))
}

#' Estimate multilocus allele frequency and COI with SNP-Slice
#'
#' Estimates multilocus allele frequencies and per-specimen COI. Requires
#' **snp.slicer** (with multi-chain sampling and `estimate`) and
#' **variantstring** 1.x (Suggests).
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele / AA-style table with counts. Default column
#'   names map AA-call fields (`aa_locus`, `aa`, `reads`). See
#'   `vignette("input-formats", package = "PGEcore")`.
#' - **`loci_groups`**: Loci-groups TSV (`group_id` plus a locus column
#'   matching `target_name_col`).
#'
#' ## Outputs
#'
#' - **`mlaf_output`**: Multilocus allele frequencies (`group_id`, `variant`,
#'   `freq`, …).
#' - **`coi_output`**: COI estimates (`specimen_name`, `coi`,
#'   `coi_cons_weighted`; uncertainty columns when `estimator = "posterior"`).
#'   `coi` counts every strain assigned to a host in the best restart.
#'   `coi_cons_weighted` pools haplotype membership across all restarts and
#'   weights each haplotype by its consensus support, which counters the
#'   dictionary over-parameterisation that inflates `coi`.
#' - **`convergence_output`**: Per-restart optimisation diagnostics
#'   (`chain_id`, `seed`, `map_logpost`, `is_best`, `map_iteration`,
#'   `final_iteration`, `plateau_frac`, `map_kstar`, `map_ktrunc`, `coi_mean`,
#'   `coi_ccc_to_best`). SNP-Slice reports the restart with the highest MAP log
#'   posterior rather than pooling chains, so between-chain R-hat and ESS do not
#'   describe its output and are not emitted.
#'
#' ## Running
#'
#' ```r
#' snpslice_wrapper(
#'   allele_table = "aa_calls.tsv",
#'   loci_groups = "loci_groups.tsv",
#'   mlaf_output = "mlaf.tsv",
#'   coi_output = "coi.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/snpslice_wrapper \
#'   --allele_table aa_calls.tsv \
#'   --loci_groups loci_groups.tsv \
#'   --mlaf_output mlaf.tsv \
#'   --coi_output coi.tsv
#' ```
#'
#' Requires **snp.slicer** and **variantstring** (Suggests).
#'
#' @param allele_table Path to allele / AA-calls TSV with counts. See *Inputs*.
#' @param loci_groups Path to loci-groups TSV. See *Inputs*.
#' @param mlaf_output Path for multilocus allele-frequency TSV. See *Outputs*.
#' @param coi_output Path for COI TSV. See *Outputs*.
#' @param convergence_output Path for per-restart optimisation-diagnostics TSV.
#'   See *Outputs*.
#' @param specimen_name_col,target_name_col,target_value_col,target_count_col
#'   Column names in `allele_table`.
#' @param loci_limit Optional cap on the number of loci.
#' @param model Observation model for SNP-Slice.
#' @param n_sample Post-burn-in MCMC iterations retained per chain.
#' @param n_burnin Burn-in iterations per chain. If `NULL`, SNP-Slice uses
#'   `floor(n_sample / 2)`.
#' @param alpha,threshold,gap SNP-Slice MCMC settings.
#' @param rho Dictionary sparsity parameter. If `NULL`, it is left to
#'   [snp.slicer::snp_slice()], which defaults to 0.5 for the categorical model
#'   and the global minor allele frequency for count models.
#' @param estimator Estimator for COI and allele frequencies: `"final_sample"`
#'   (default, matching the SNP-Slice paper), `"map"`, or `"posterior"` (the
#'   posterior mean, which also yields uncertainty columns).
#' @param n_chains Number of independent MCMC chains. More than one is required
#'   for the Gelman-Rubin R-hat diagnostic.
#' @param threads Cores used to run chains simultaneously (capped at `n_chains`).
#' @param verbose Verbose SNP-Slice output.
#' @param seed Random seed.
#'
#' @return Invisibly, a list with `mlaf`, `coi`, and `convergence` tibbles.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
snpslice_wrapper <- function(allele_table,
                             loci_groups,
                             mlaf_output,
                             coi_output,
                             convergence_output = "convergence_diag.tsv",
                             specimen_name_col = "specimen_name",
                             target_name_col = "aa_locus",
                             target_value_col = "aa",
                             target_count_col = "reads",
                             loci_limit = NULL,
                             model = "negative_binomial",
                             n_sample = 10000L,
                             n_burnin = NULL,
                             alpha = 2.6,
                             rho = NULL,
                             threshold = 0.001,
                             gap = NULL,
                             estimator = "final_sample",
                             n_chains = 3L,
                             threads = 1L,
                             verbose = FALSE,
                             seed = 1L) {
  check_suggested_pkg("snp.slicer", "SNP-Slice via snpslice_wrapper()")
  check_variantstring_v1("variant strings via snpslice_wrapper()")

  required <- list(
    allele_table = allele_table,
    loci_groups = loci_groups,
    mlaf_output = mlaf_output,
    coi_output = coi_output
  )
  missing <- names(required)[vapply(required, is.null, logical(1))]
  if (length(missing) > 0) {
    stop(
      "missing the following arguments: ",
      paste(paste0("--", missing), collapse = ", "),
      call. = FALSE
    )
  }
  valid_estimators <- c("final_sample", "map", "posterior")
  if (!estimator %in% valid_estimators) {
    stop(
      "--estimator must be one of: ",
      paste(valid_estimators, collapse = ", "),
      call. = FALSE
    )
  }

  set.seed(seed)
  allele_tbl <- create_snpslice_allele_table_input(
    allele_table,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col,
    target_count_col = target_count_col
  )
  loci_groups <- create_snpslice_loci_group_input(
    loci_groups,
    allele_tbl,
    target_name_col = target_name_col
  )

  if (!is.null(loci_limit)) {
    if (dplyr::n_distinct(allele_tbl$target_name) > loci_limit) {
      loci_oi <- unique(unlist(loci_groups))
      n_loci_select <- max(0L, loci_limit - length(loci_oi))
      if (n_loci_select == 0L) {
        message(
          "Note: loci_limit (", loci_limit, ") is <= the number of group loci (",
          length(loci_oi), "); no extra loci will be sampled."
        )
      }
      # SNP-Slice keeps only targets with at most two alleles, so the extra
      # loci are drawn from the biallelic targets alone. Monomorphic targets
      # are excluded because they carry no allelic variation.
      biallelic_trgs <- allele_tbl |>
        dplyr::group_by(.data$target_name) |>
        dplyr::filter(dplyr::n_distinct(.data$target_value) == 2) |>
        dplyr::pull("target_name") |>
        unique()
      candidates <- setdiff(biallelic_trgs, loci_oi)
      if (n_loci_select > length(candidates)) {
        message(
          "Note: only ", length(candidates), " biallelic target(s) are ",
          "available to sample; requested ", n_loci_select, "."
        )
        n_loci_select <- length(candidates)
      }
      loci_random <- sample(candidates, n_loci_select)
      loci_selected <- c(loci_oi, loci_random)
      allele_tbl <- dplyr::filter(allele_tbl, .data$target_name %in% loci_selected)
    }
  }

  snpslice_args <- list(
    allele_tbl,
    model = model,
    n_sample = n_sample,
    n_burnin = n_burnin,
    alpha = alpha,
    threshold = threshold,
    gap = gap,
    n_chains = n_chains,
    n_cores = threads,  # snp.slicer arg name
    seed = seed,
    # Retained per-iteration samples are only needed by the "posterior"
    # estimator.
    store_mcmc = identical(estimator, "posterior"),
    verbose = verbose,
    specimen_id_col = "specimen_name",
    target_id_col = "target_name",
    target_value_col = "target_value",
    target_count_col = "target_count"
  )
  # rho is only passed when supplied. Passing rho = NULL explicitly would
  # bypass the model-specific default in snp_slice() (0.5 for categorical, the
  # minor allele frequency for count models) and force the minor allele
  # frequency for every model, which is wrong for categorical data.
  if (!is.null(rho)) {
    snpslice_args$rho <- rho
  }
  snpslice_res <- do.call(snp.slicer::snp_slice, snpslice_args)

  mlaf <- prepare_snpslice_af_output(snpslice_res, loci_groups, estimator)
  readr::write_tsv(mlaf, mlaf_output)
  coi <- prepare_snpslice_coi_output(snpslice_res, specimen_name_col, estimator)
  readr::write_tsv(coi, coi_output)
  convergence <- prepare_snpslice_optim_output(snpslice_res)
  readr::write_tsv(convergence, convergence_output)
  invisible(list(mlaf = mlaf, coi = coi, convergence = convergence))
}
