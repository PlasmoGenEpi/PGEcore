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
#' @param allow_multiallelic Keep targets with more than two alleles. The
#'   biallelic observation models drop them (SNP-Slice cannot represent a
#'   third allele), so by default they are excluded from every group with a
#'   warning; the multinomial model keeps them.
#' @keywords internal
create_snpslice_loci_group_input <- function(loci_groups_path,
                                             allele_table,
                                             target_name_col = "target_name",
                                             allow_multiallelic = FALSE) {
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
    non_biallelic_trgs <- if (allow_multiallelic) {
      character()
    } else {
      allele_table |>
        dplyr::filter(.data$target_name %in% loci_groups[[lg]]) |>
        dplyr::group_by(.data$target_name) |>
        dplyr::filter(dplyr::n_distinct(.data$target_value) > 2) |>
        dplyr::pull("target_name") |>
        unique()
    }
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

#' Collapse identical dictionary rows of a SNP-Slice allocation
#'
#' Ju et al. (2024) count identical SNP haplotypes once when computing MOI,
#' allele frequencies and heterozygosity: duplicated dictionary rows are
#' removed before any estimate is formed. SNP-Slice can hold two active
#' strains with the same haplotype, and a specimen assigned both would then
#' contribute two to its COI and two strain copies to every frequency. This
#' merges such strains: a specimen carries the merged strain if it carried
#' any of the duplicates.
#'
#' @param A Allocation matrix, specimens x strains.
#' @param D Dictionary matrix, strains x targets.
#' @return List with the merged `A` and `D`; strains keep first-appearance
#'   order.
#' @keywords internal
snpslice_dedup_matrices <- function(A, D) {
  if (nrow(D) == 0L) {
    return(list(A = A, D = D))
  }
  key <- apply(D, 1L, paste0, collapse = ",")
  first <- match(unique(key), key)
  members <- split(seq_along(key), factor(key, levels = key[first]))
  A_merged <- vapply(members, function(k) {
    as.numeric(rowSums(A[, k, drop = FALSE]) > 0)
  }, numeric(nrow(A)))
  A_merged <- matrix(A_merged, nrow = nrow(A), dimnames = list(rownames(A), NULL))
  list(A = A_merged, D = D[first, , drop = FALSE])
}

#' Apply [snpslice_dedup_matrices()] to every estimate a chain carries
#'
#' @param chain A single-chain `snp_slice_results` object from
#'   [snp.slicer::get_chain()].
#' @return The chain with its MAP, final-sample and (if present) stored MCMC
#'   sample matrices de-duplicated.
#' @keywords internal
snpslice_dedup_chain <- function(chain) {
  m <- snpslice_dedup_matrices(chain$map_allocation_matrix, chain$map_dictionary_matrix)
  chain$map_allocation_matrix <- m$A
  chain$map_dictionary_matrix <- m$D
  if (!is.null(chain$final_allocation_matrix)) {
    f <- snpslice_dedup_matrices(chain$final_allocation_matrix, chain$final_dictionary_matrix)
    chain$final_allocation_matrix <- f$A
    chain$final_dictionary_matrix <- f$D
  }
  if (!is.null(chain$mcmc_samples)) {
    chain$mcmc_samples <- lapply(chain$mcmc_samples, function(s) {
      d <- snpslice_dedup_matrices(s$A, s$D)
      s$A <- d$A
      s$D <- d$D
      s
    })
  }
  chain
}

#' Split a SNP-Slice result into single-chain objects
#'
#' @param snpslice_res Result of [snp.slicer::snp_slice()].
#' @param dedup_haplotypes Collapse identical haplotypes in every chain; see
#'   [snpslice_dedup_matrices()].
#' @return List with `chains` (single-chain results in chain order) and
#'   `best`, the index of the chain SNP-Slice reports (highest MAP log
#'   posterior).
#' @keywords internal
snpslice_split_chains <- function(snpslice_res, dedup_haplotypes = FALSE) {
  n <- if (is.null(snpslice_res$chains)) 1L else length(snpslice_res$chains)
  chains <- lapply(seq_len(n), function(i) snp.slicer::get_chain(snpslice_res, i))
  if (dedup_haplotypes) {
    chains <- lapply(chains, snpslice_dedup_chain)
  }
  best <- if (is.null(snpslice_res$best_chain)) 1L else as.integer(snpslice_res$best_chain)
  list(chains = chains, best = best)
}

#' Format SNP-Slice allele frequencies with variantstring names
#'
#' `freq` is the estimate from the chain SNP-Slice reports (highest MAP log
#' posterior). `averaged_freq` follows Ju et al. (2024), who report the mean
#' of the per-chain estimates over independently initialised chains: it is
#' the mean over all chains of that chain's frequency for the haplotype, with
#' a haplotype a chain never produced counting as zero there. A row is kept
#' when either column is positive.
#'
#' @param chains,best Output of [snpslice_split_chains()].
#' @keywords internal
prepare_snpslice_af_output <- function(chains, best, loci_groups, estimator) {
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
      dplyr::filter(.data$frequency > 0 | .data$averaged_freq > 0) |>
      dplyr::mutate(
        allele = lapply(
          .data$allele,
          prep_variantstring_input,
          loci_groups[[group_id]]
        )
      ) |>
      dplyr::mutate(allele = variantstring::long_to_variant(.data$allele))
  }

  # One frequency table per group per chain. Every chain enumerates the same
  # haplotype rows for a group, so they can be joined on `allele`.
  per_chain <- lapply(chains, function(ch) {
    snp.slicer::calculate_allele_frequencies_by_sets(ch, loci_groups, estimate = estimator)
  })
  by_group <- lapply(names(loci_groups), function(g) {
    tab <- tibble::as_tibble(per_chain[[best]][[g]])
    freq_mat <- vapply(per_chain, function(pc) {
      pc[[g]]$frequency[match(tab$allele, pc[[g]]$allele)]
    }, numeric(nrow(tab)))
    freq_mat <- matrix(freq_mat, nrow = nrow(tab))
    freq_mat[is.na(freq_mat)] <- 0
    tab$averaged_freq <- rowMeans(freq_mat)
    tab
  })
  names(by_group) <- names(loci_groups)

  af_tib <- tibble::tibble(
    group_id = names(by_group),
    af_tib = unname(by_group)
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
    dplyr::rename(variant = "allele", freq = "frequency") |>
    dplyr::relocate("averaged_freq", .after = "freq")

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
#' Pools the strain assignments from several independent SNP-Slice restarts
#' into one per-host complexity of infection (COI) that discounts strains
#' the restarts do not agree on.
#'
#' @section Why this exists:
#' SNP-Slice fits allele frequencies well but over-parameterizes the strain
#' dictionary to do so: it adds many low-support strains, often carried by a
#' single host. Each such strain adds a full +1 to that host's COI while
#' contributing almost nothing to the frequencies, so the raw row sum of the
#' allocation matrix over-counts COI. On top of that, restarts are independent
#' optimizations that land on different dictionaries, so any single restart's
#' assignments are partly noise. This function addresses both problems by
#' asking, per host and per strain, how consistently that assignment is
#' recovered across the different chains, and weighting the count accordingly.
#'
#' @section How it works:
#' 1. **Match strains across restarts by sequence.** A strain index means
#'    nothing across restarts, so haplotypes are keyed by their dictionary
#'    row (the concatenated allele string). Duplicate rows within a restart
#'    are collapsed before counting.
#' 2. **Build a membership matrix.** `membership[i, h]` is the fraction of
#'    restarts in which host `i` carries haplotype `h`. Assignments
#'    recovered by every restart score 1; those seen in one restart out of
#'    ten score 0.1.
#' 3. **Weight each haplotype by cohort support.** Support is
#'    `colSums(membership)`, the expected number of hosts carrying `h`.
#'    Each haplotype gets weight `1 - exp(-support / mean(support))`.
#'    This is a smooth discount rather than a threshold: a strain with no
#'    support contributes 0, one of average commonness contributes about
#'    0.63, and a common strain contributes fully. Hard thresholds were
#'    tested and rejected because the best cutoff differed by population
#'    and a small wobble in support flipped whole-strain counts.
#' 4. **Sum and floor.** Host COI is `membership %*% weights`, floored at 1
#'    so every host is counted as at least one infection.
#'
#' Normalizing support by its mean, rather than a fixed host count, is what
#' keeps the estimate from drifting with panel or cohort size. More loci
#' resolve more strains, which lowers the mean support, so the same absolute
#' support earns a higher weight. More hosts raise the mean, so the same
#' support earns a lower weight. Both adjustments are the desired direction.
#'
#' In benchmarking on simulated populations this consensus estimate beat
#' the same weighting applied to a single restart, nearly removed the loci
#' drift, and was more reproducible between independent seeds. The gain
#' saturates at roughly three restarts.
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
#' Three per-specimen columns:
#' - `coi`: strains assigned in the reported chain (highest MAP log
#'   posterior), the row sum of its allocation matrix.
#' - `coi_cons_weighted`: [snpslice_consensus_coi()] across chains.
#' - `coi_chain_mean`: the mean over all chains of that chain's `coi`, which
#'   is how Ju et al. (2024) report MOI (the average of the final-sample
#'   estimate over independently initialised chains).
#'
#' @param chains,best Output of [snpslice_split_chains()].
#' @keywords internal
prepare_snpslice_coi_output <- function(chains,
                                        best,
                                        specimen_name_col,
                                        estimator) {
  coi_tib <- snp.slicer::calculate_individual_coi(
    chains[[best]],
    estimate = estimator
  ) |>
    dplyr::select(-"host_index") |>
    dplyr::rename(!!specimen_name_col := "host_id", coi = "coi_estimate")
  if (!identical(estimator, "posterior")) {
    coi_tib <- dplyr::select(coi_tib, -"coi_sd", -"coi_lower", -"coi_upper")
  }
  point_estimate <- if (identical(estimator, "posterior")) "map" else estimator
  consensus <- snpslice_consensus_coi(chains, point_estimate)
  if (length(consensus) != nrow(coi_tib)) {
    stop("Consensus COI length does not match the per-host COI table.",
         call. = FALSE)
  }
  per_chain_coi <- vapply(chains, function(ch) {
    snp.slicer::calculate_individual_coi(ch, estimate = estimator)$coi_estimate
  }, numeric(nrow(coi_tib)))
  coi_tib$coi_cons_weighted <- consensus
  coi_tib$coi_chain_mean <- rowMeans(matrix(per_chain_coi, nrow = nrow(coi_tib)))
  dplyr::relocate(coi_tib, "coi_cons_weighted", "coi_chain_mean", .after = "coi")
}

#' Lin's concordance correlation coefficient
#'
#' Computes Lin's (1989) concordance correlation coefficient (CCC) for
#' agreement between two sets of paired measurements. Unlike Pearson's
#' correlation, the CCC combines precision (tightness of the points about
#' their best-fit line) and accuracy (how far that line deviates from the
#' 45-degree line of perfect concordance), so it measures reproducibility
#' rather than linear association. Values range from -1 to 1, with 1
#' indicating perfect agreement. Used here to compare per-host COI estimates
#' between SNP-Slice restarts.
#'
#' Only pairs where both `x` and `y` are non-missing are used. Returns
#' `NA_real_` for fewer than three complete pairs and 1 when both vectors are
#' constant and equal.
#'
#' @param x Numeric vector, the first set of measurements.
#' @param y Numeric vector, the second set of measurements, same length as
#'   `x`.
#'
#' @return A single numeric value, the concordance correlation coefficient.
#'
#' @references
#' Lin L (1989). A concordance correlation coefficient to evaluate
#' reproducibility. *Biometrics* 45: 255-268.
#'
#' Lin L (2000). A note on the concordance correlation coefficient.
#' *Biometrics* 56: 324-325.
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

#' Format SNP-Slice per-restart optimization diagnostics
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
#' **snp.slicer** (with multi-chain sampling and `estimate`) and, when loci
#' groups are supplied, **variantstring** 1.x (Suggests).
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele / AA-style table with counts. Default column
#'   names map AA-call fields (`aa_locus`, `aa`, `reads`). See
#'   `vignette("input-formats", package = "PGEcore")`.
#' - **`loci_groups`** (optional): Loci-groups TSV (`group_id` plus a locus
#'   column matching `target_name_col`). Supply together with `mlaf_output`.
#'   Omit both to run in COI-only mode.
#'
#' ## Multi-allelic targets
#'
#' The biallelic observation models (`"negative_binomial"`, `"binomial"`,
#' `"poisson"`, `"categorical"`) cannot represent a third allele, so any target
#' with more than two alleles is dropped from its group with a warning, and
#' extra loci drawn under `loci_limit` come from biallelic targets only. With
#' `model = "multinomial"` every allele at every target is kept, groups with
#' multi-allelic codons (for example dhfr 108 or dhps 540) are estimated in
#' full, and extra loci are drawn from all polymorphic targets. The
#' multinomial dictionary prior is set by `dict_prior`; `rho` is ignored.
#'
#' ## Outputs
#'
#' - **`mlaf_output`** (only with `loci_groups`): Multilocus allele frequencies
#'   (`group_id`, `variant`, `freq`, `averaged_freq`, …). `freq` comes from the
#'   chain SNP-Slice reports (highest MAP log posterior); `averaged_freq` is the
#'   mean of the per-chain estimates over all chains, the estimator of
#'   Ju et al. (2024).
#' - **`coi_output`**: COI estimates (`specimen_name`, `coi`,
#'   `coi_cons_weighted`, `coi_chain_mean`; uncertainty columns when
#'   `estimator = "posterior"`).
#'   `coi` counts every strain assigned to a host in the best restart.
#'   `coi_cons_weighted` pools haplotype membership across all restarts and
#'   weights each haplotype by its consensus support, which counters the
#'   dictionary over-parameterisation that inflates `coi`. `coi_chain_mean`
#'   is the mean of `coi` over all restarts, as Ju et al. (2024) report MOI.
#'
#' With `dedup_haplotypes = TRUE` (the default), strains with identical
#' haplotypes are merged in every chain before any of these are computed, so
#' a specimen assigned two copies of the same haplotype counts it once, as in
#' Ju et al. (2024). `FALSE` restores the earlier raw counts.
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
#'   --allele_table ./inst/extdata/example_aa_calls.tsv \
#'   --loci_groups loci_groups.tsv \
#'   --mlaf_output mlaf.tsv \
#'   --coi_output coi.tsv
#' ```
#'
#' COI only:
#'
#' ```bash
#' Rscript exec/snpslice_wrapper \
#'   --allele_table aa_calls.tsv \
#'   --loci_limit 100 \
#'   --coi_output coi.tsv
#' ```
#'
#' Requires **snp.slicer**, plus **variantstring** when `loci_groups` is
#' supplied (Suggests).
#'
#' @param allele_table Path to allele / AA-calls TSV with counts. See *Inputs*.
#' @param loci_groups Optional path to loci-groups TSV. Must be supplied
#'   together with `mlaf_output`. See *Inputs*.
#' @param mlaf_output Optional path for multilocus allele-frequency TSV. Must
#'   be supplied together with `loci_groups`. See *Outputs*.
#' @param coi_output Path for COI TSV. See *Outputs*.
#' @param convergence_output Path for per-restart optimization-diagnostics TSV.
#'   See *Outputs*.
#' @param specimen_name_col,target_name_col,target_value_col,target_count_col
#'   Column names in `allele_table`.
#' @param loci_limit Optional cap on the number of loci. With `loci_groups`,
#'   the group loci are always kept and the remainder is filled with random
#'   loci; without, every selected locus is random. Random loci are drawn from
#'   the biallelic targets, or from any polymorphic target when `model` is
#'   `"multinomial"`. If `NULL`, all loci are used.
#' @param model Observation model for SNP-Slice: `"negative_binomial"`
#'   (default), `"binomial"`, `"poisson"`, `"categorical"`, or
#'   `"multinomial"`. See *Multi-allelic targets*.
#' @param dict_prior Dictionary prior for the multinomial model, passed to
#'   [snp.slicer::snp_slice()]: `"empirical"` (default; pooled allele read
#'   fractions with a pseudocount) or `"uniform"`. Ignored by other models.
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
#' @param dedup_haplotypes Merge strains with identical haplotypes in every
#'   chain before computing COI and frequencies. See *Outputs*.
#' @param n_chains Number of independent MCMC chains. More than one is required
#'   for the Gelman-Rubin R-hat diagnostic.
#' @param threads Cores used to run chains simultaneously (capped at `n_chains`).
#' @param verbose Verbose SNP-Slice output.
#' @param seed Random seed.
#'
#' @return Invisibly, a list with `mlaf`, `coi`, and `convergence` tibbles.
#'   `mlaf` is `NULL` when `loci_groups` is not supplied.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
snpslice_wrapper <- function(allele_table,
                             coi_output,
                             loci_groups = NULL,
                             mlaf_output = NULL,
                             convergence_output = "convergence_diag.tsv",
                             specimen_name_col = "specimen_name",
                             target_name_col = "aa_locus",
                             target_value_col = "aa",
                             target_count_col = "reads",
                             loci_limit = NULL,
                             model = "negative_binomial",
                             dict_prior = "empirical",
                             n_sample = 10000L,
                             n_burnin = NULL,
                             alpha = 2.6,
                             rho = NULL,
                             threshold = 0.001,
                             gap = NULL,
                             estimator = "final_sample",
                             dedup_haplotypes = TRUE,
                             n_chains = 3L,
                             threads = 1L,
                             verbose = FALSE,
                             seed = 1L) {
  check_suggested_pkg("snp.slicer", "SNP-Slice via snpslice_wrapper()")
  if (!is.null(loci_groups)) {
    check_variantstring_v1("variant strings via snpslice_wrapper()")
  }

  required <- list(
    allele_table = allele_table,
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
  if (is.null(loci_groups) != is.null(mlaf_output)) {
    stop(
      "--loci_groups and --mlaf_output must be supplied together. Omit both ",
      "to estimate COI only.",
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

  valid_models <- c(
    "negative_binomial", "binomial", "poisson", "categorical", "multinomial"
  )
  if (!model %in% valid_models) {
    stop(
      "--model must be one of: ",
      paste(valid_models, collapse = ", "),
      call. = FALSE
    )
  }
  # Only the multinomial model can carry more than two alleles per target.
  multiallelic <- identical(model, "multinomial")

  set.seed(seed)
  allele_tbl <- create_snpslice_allele_table_input(
    allele_table,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col,
    target_count_col = target_count_col
  )
  if (!is.null(loci_groups)) {
    loci_groups <- create_snpslice_loci_group_input(
      loci_groups,
      allele_tbl,
      target_name_col = target_name_col,
      allow_multiallelic = multiallelic
    )
  }

  if (!is.null(loci_limit)) {
    if (dplyr::n_distinct(allele_tbl$target_name) > loci_limit) {
      # Without loci groups (COI-only mode) loci_oi is empty, so every
      # selected locus is drawn at random.
      loci_oi <- unique(unlist(loci_groups))
      n_loci_select <- max(0L, loci_limit - length(loci_oi))
      if (n_loci_select == 0L) {
        message(
          "Note: loci_limit (", loci_limit, ") is <= the number of group loci (",
          length(loci_oi), "); no extra loci will be sampled."
        )
      }
      # The biallelic models keep only targets with at most two alleles, so
      # their extra loci are drawn from the biallelic targets alone; the
      # multinomial model can draw from any polymorphic target. Monomorphic
      # targets are excluded because they carry no allelic variation.
      n_alleles_by_target <- allele_tbl |>
        dplyr::group_by(.data$target_name) |>
        dplyr::summarise(
          n_alleles = dplyr::n_distinct(.data$target_value),
          .groups = "drop"
        )
      eligible <- if (multiallelic) {
        n_alleles_by_target$n_alleles >= 2
      } else {
        n_alleles_by_target$n_alleles == 2
      }
      candidate_kind <- if (multiallelic) "polymorphic" else "biallelic"
      candidates <- setdiff(n_alleles_by_target$target_name[eligible], loci_oi)
      if (n_loci_select > length(candidates)) {
        message(
          "Note: only ", length(candidates), " ", candidate_kind,
          " target(s) are available to sample; requested ", n_loci_select, "."
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
  # dict_prior is a multinomial-only argument; the other models' loaders
  # would reject it.
  if (multiallelic) {
    snpslice_args$dict_prior <- dict_prior
  }
  snpslice_res <- do.call(snp.slicer::snp_slice, snpslice_args)

  split <- snpslice_split_chains(snpslice_res, dedup_haplotypes = dedup_haplotypes)
  mlaf <- NULL
  if (!is.null(loci_groups)) {
    mlaf <- prepare_snpslice_af_output(split$chains, split$best, loci_groups, estimator)
    readr::write_tsv(mlaf, mlaf_output)
  }
  coi <- prepare_snpslice_coi_output(split$chains, split$best, specimen_name_col, estimator)
  readr::write_tsv(coi, coi_output)
  convergence <- prepare_snpslice_optim_output(snpslice_res)
  readr::write_tsv(convergence, convergence_output)
  invisible(list(mlaf = mlaf, coi = coi, convergence = convergence))
}
