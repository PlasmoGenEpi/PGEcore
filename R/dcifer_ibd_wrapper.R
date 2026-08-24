#' Read allele frequencies into Dcifer list format
#'
#' @keywords internal
create_dcifer_allele_freq_input <- function(allele_freq_path,
                                            allele_list,
                                            target_name_col = "target_name",
                                            target_value_col = "seq") {
  allele_freqs <- readr::read_tsv(
    allele_freq_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      freq = readr::col_double()
    ),
    progress = FALSE
  ) |>
    dplyr::select(dplyr::all_of(c(target_name_col, target_value_col)), "freq") |>
    dplyr::rename(
      target_name = dplyr::all_of(target_name_col),
      target_value = dplyr::all_of(target_value_col)
    )

  rules <- validate::validator(
    is.character(target_name),
    is.character(target_value),
    is.double(freq),
    !is.na(target_name),
    !is.na(target_value),
    !is.na(freq)
  )
  stop_on_validate_fails(allele_freqs, rules, "allele_freq_table")

  allele_freq_alleles <- allele_freqs |>
    tidyr::unite("allele", "target_name", "target_value", sep = ":") |>
    dplyr::pull("allele")

  collapsed <- purrr::list_c(allele_list)
  allele_table_alleles <- tibble::tibble(
    target_name = names(collapsed),
    target_values = unname(collapsed)
  ) |>
    dplyr::mutate(target_values = lapply(.data$target_values, names)) |>
    tidyr::unnest("target_values") |>
    tidyr::unite("alleles", "target_name", "target_values", sep = ":") |>
    dplyr::pull("alleles")

  missing_in_freq <- setdiff(allele_table_alleles, allele_freq_alleles)
  if (length(missing_in_freq) > 0) {
    stop(
      "The following alleles are in the allele table and not in the provided ",
      "allele frequencies: ",
      paste(missing_in_freq, collapse = " "),
      call. = FALSE
    )
  }
  extra_in_freq <- setdiff(allele_freq_alleles, allele_table_alleles)
  if (length(extra_in_freq) > 0) {
    stop(
      "The following alleles are in the provided allele frequencies and not ",
      "in the allele table: ",
      paste(extra_in_freq, collapse = " "),
      call. = FALSE
    )
  }

  dcifer::formatAfreq(
    allele_freqs,
    lvar = "target_name",
    avar = "target_value",
    fvar = "freq"
  )
}

#' Read specimen metadata for population-specific Dcifer runs
#'
#' @keywords internal
create_dcifer_specimen_metadata_input <- function(specimen_metadata_path,
                                                  coi,
                                                  specimen_name_col = "specimen_name",
                                                  pop_name_col = "population") {
  specimen_metadata <- readr::read_tsv(
    specimen_metadata_path,
    col_types = do.call(
      readr::cols,
      c(
        stats::setNames(list(readr::col_character()), specimen_name_col),
        list(.default = readr::col_character())
      )
    ),
    progress = FALSE
  ) |>
    dplyr::select(dplyr::all_of(c(specimen_name_col, pop_name_col))) |>
    dplyr::rename(
      specimen_name = dplyr::all_of(specimen_name_col),
      population = dplyr::all_of(pop_name_col)
    )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(population),
    !is.na(specimen_name),
    !is.na(population)
  )
  stop_on_validate_fails(specimen_metadata, rules, "specimen_metadata")

  specimens_notinmeta <- setdiff(names(coi), specimen_metadata$specimen_name)
  if (length(specimens_notinmeta) > 0) {
    stop(
      "Some specimen IDs are missing from metadata: ",
      paste(specimens_notinmeta, collapse = ","),
      call. = FALSE
    )
  }
  specimen_metadata
}

#' Mean of two Dcifer allele-frequency lists
#'
#' @keywords internal
dcifer_allele_freq_mean <- function(af1, af2) {
  mean_af <- list()
  loci <- unique(c(names(af1), names(af2)))
  for (l in loci) {
    alleles <- unique(c(names(af1[[l]]), names(af2[[l]])))
    locus_mean_freqs <- numeric()
    for (a in alleles) {
      af1_af2 <- c(af1[[l]][a], af2[[l]][a])
      af1_af2[is.na(af1_af2)] <- 0
      locus_mean_freqs[a] <- mean(af1_af2)
    }
    mean_af[[l]] <- locus_mean_freqs
  }
  mean_af
}

#' Run ibdPair in parallel
#'
#' @keywords internal
run_ibdpair <- function(dsmp,
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
  check_suggested_pkgs(
    c("dcifer", "foreach", "iterators", "doParallel", "parallelly"),
    "IBD relatedness via run_ibdpair()"
  )
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  if (confint) {
    mnewton <- FALSE
    tol <- NULL
  } else {
    mnewton <- TRUE
    tol <- 1 / nr
  }
  if (!mnewton) {
    if (!inherits(reval, "matrix")) {
      reval <- dcifer::generateReval(M = 1, rval = reval, nr = nr)
    }
    neval <- ncol(reval)
    logr <- dcifer::logReval(reval, M = 1, neval = neval)
  } else {
    neval <- NULL
    logr <- NULL
  }
  inull <- if (mnewton || !pval) {
    NULL
  } else {
    which.min(abs(reval - rnull))
  }
  afreq <- lapply(afreq, log)
  nloc <- length(afreq)

  if (identical(rnull, 0L)) {
    side <- "right"
  } else if (identical(rnull, 1L)) {
    side <- "left"
  } else {
    side <- "two-sided"
  }

  if (is.null(sample_pairs)) {
    sample_pairs_matrix <- utils::combn(names(dsmp), 2)
    sample_pairs <- tibble::tibble(
      sample_a = sample_pairs_matrix[1, ],
      sample_b = sample_pairs_matrix[2, ]
    )
  }

  if (is.null(total_cores)) {
    total_cores <- parallelly::availableCores() - 1
  }

  if (is.null(parallel::getDefaultCluster())) {
    if (verbose) {
      cl <- parallel::makeCluster(total_cores, outfile = "")
    } else {
      cl <- parallel::makeCluster(total_cores)
    }
    parallel::setDefaultCluster(cl)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- parallel::getDefaultCluster()
  }

  res <- foreach::foreach(
    i = 1:total_cores,
    .combine = rbind,
    .packages = c("dcifer", "foreach", "iterators")
  ) %dopar% {
    total_pairs <- nrow(sample_pairs)
    begin_idx <- floor(((i - 1) * total_pairs / total_cores) + 1)
    end_idx <- floor((i * total_pairs / total_cores))
    pairs <- sample_pairs[begin_idx:end_idx, ]
    foreach::foreach(
      pair = iterators::iter(pairs, by = "row"),
      .combine = rbind,
      .verbose = verbose
    ) %do% {
      sample_a <- pair$sample_a
      sample_b <- pair$sample_b
      rxy <- dcifer::ibdPair(
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
      tibble::tibble(
        sample_a = sample_a,
        sample_b = sample_b,
        estimate = rxy$rhat,
        p_value = rxy$pval,
        CI_lower = range(rxy$confreg)[1],
        CI_upper = range(rxy$confreg)[2]
      )
    }
  }
  parallel::stopCluster(cl)
  parallel::setDefaultCluster(NULL)
  res
}

#' Run ibdEstM in parallel
#'
#' @keywords internal
run_ibdestm <- function(dsmp,
                        coi,
                        afreq,
                        sample_pairs = NULL,
                        pval = TRUE,
                        confint = FALSE,
                        rnull = 0,
                        alpha = 0.05,
                        total_cores = NULL,
                        verbose = FALSE) {
  check_suggested_pkgs(
    c("dcifer", "foreach", "iterators", "doParallel", "parallelly"),
    "IBD relatedness via run_ibdestm()"
  )
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  nrs <- c(1e3, 1e2, 32, 16, 12, 10)
  revals <- mapply(dcifer::generateReval, 1:6, nr = nrs)
  afreq <- lapply(afreq, log)
  nloc <- length(afreq)
  if (is.null(sample_pairs)) {
    sample_pairs_matrix <- utils::combn(names(dsmp), 2)
    sample_pairs <- tibble::tibble(
      sample_a = sample_pairs_matrix[1, ],
      sample_b = sample_pairs_matrix[2, ]
    )
  }

  if (identical(rnull, 0L)) {
    side <- "right"
  } else if (identical(rnull, 1L)) {
    side <- "left"
  } else {
    side <- "two-sided"
  }

  if (is.null(total_cores)) {
    total_cores <- parallelly::availableCores() - 1
  }
  if (is.null(parallel::getDefaultCluster())) {
    if (verbose) {
      cl <- parallel::makeCluster(total_cores, outfile = "")
    } else {
      cl <- parallel::makeCluster(total_cores)
    }
    parallel::setDefaultCluster(cl)
    doParallel::registerDoParallel(cl)
  } else {
    cl <- parallel::getDefaultCluster()
  }

  res <- foreach::foreach(
    i = 1:total_cores,
    .combine = rbind,
    .packages = c("dcifer", "foreach", "iterators")
  ) %dopar% {
    total_pairs <- nrow(sample_pairs)
    begin_idx <- floor(((i - 1) * total_pairs / total_cores) + 1)
    end_idx <- floor((i * total_pairs / total_cores))
    pairs <- sample_pairs[begin_idx:end_idx, ]
    foreach::foreach(
      pair = iterators::iter(pairs, by = "row"),
      .combine = rbind,
      .verbose = verbose
    ) %do% {
      sample_a <- pair$sample_a
      sample_b <- pair$sample_b
      rxy <- dcifer::ibdEstM(
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
      tibble::tibble(
        sample_a = sample_a,
        sample_b = sample_b,
        estimate = estimate,
        strain_pair = seq_along(estimate),
        p_value = rxy$pval,
        CI_lower = range(rxy$confreg)[1],
        CI_upper = range(rxy$confreg)[2]
      )
    }
  }
  parallel::stopCluster(cl)
  parallel::setDefaultCluster(NULL)
  res
}

#' Run Dcifer with population-specific allele frequencies
#'
#' @keywords internal
run_dcifer_bypop <- function(allele_table_path,
                             coi,
                             specimen_metadata,
                             specimen_name_col = "specimen_name",
                             target_name_col = "target_name",
                             target_value_col = "seq",
                             use_estm = FALSE,
                             ...) {
  allele_table <- create_dcifer_allele_table_input(
    allele_table_path,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col
  )
  alleles_w_specimen_meta <- allele_table |>
    dplyr::left_join(specimen_metadata, by = "specimen_name")

  pops <- unique(alleles_w_specimen_meta$population)
  allele_freq_lists <- list()
  dcifer_res <- list()
  for (pop_oi in pops) {
    alleles_filtered <- alleles_w_specimen_meta |>
      dplyr::filter(.data$population == pop_oi)
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
    if (use_estm) {
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

  if (length(pops) > 1) {
    pop_combos <- utils::combn(pops, 2)
    for (cb in seq_len(ncol(pop_combos))) {
      pop_combo <- pop_combos[, cb]
      pop_oi <- paste(sort(pop_combo), collapse = "_")
      pop_combo_af <- dcifer_allele_freq_mean(
        allele_freq_lists[[pop_combo[1]]],
        allele_freq_lists[[pop_combo[2]]]
      )
      matched <- alleles_w_specimen_meta |>
        dplyr::filter(.data$population %in% pop_combo) |>
        dcifer::formatDat(
          svar = "specimen_name",
          lvar = "target_name",
          avar = "target_value"
        ) |>
        dcifer::matchAfreq(pop_combo_af)
      pop_combo_alleles <- matched$dsmp
      pop_combo_coi <- coi[names(pop_combo_alleles)]
      samples_pop_a <- unique(
        alleles_w_specimen_meta$specimen_name[
          alleles_w_specimen_meta$population == pop_combo[1]
        ]
      )
      samples_pop_b <- unique(
        alleles_w_specimen_meta$specimen_name[
          alleles_w_specimen_meta$population == pop_combo[2]
        ]
      )
      sample_pairs <- expand.grid(samples_pop_a, samples_pop_b) |>
        tibble::as_tibble() |>
        dplyr::rename(sample_a = "Var1", sample_b = "Var2") |>
        dplyr::mutate(
          sample_a = as.character(.data$sample_a),
          sample_b = as.character(.data$sample_b)
        )
      if (use_estm) {
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

  dplyr::bind_rows(dcifer_res)
}

#' Write Dcifer IBD results
#'
#' @keywords internal
write_dcifer_ibd_output <- function(dcifer_results, out_path) {
  dcifer_results |>
    dplyr::rename(
      specimen_name_a = "sample_a",
      specimen_name_b = "sample_b",
      btwn_host_rel = "estimate"
    ) |>
    readr::write_tsv(out_path)
}

#' Estimate IBD-based relatedness with Dcifer
#'
#' File-oriented entry point used by the `dcifer_ibd_wrapper` CLI. Optional
#' **dcifer** plus parallel helpers (**doParallel**, **parallelly**, **foreach**,
#' **iterators**) must be installed separately.
#'
#' @param allele_table Path to allele TSV.
#' @param btwn_host_rel_output Path for relatedness TSV.
#' @param coi_table Optional COI TSV.
#' @param allele_freq_table Optional allele-frequency TSV.
#' @param specimen_name_col,target_name_col,target_value_col Column names.
#' @param specimen_metadata Optional metadata TSV.
#' @param pop_name_col Optional population column in metadata.
#' @param rnull Relatedness null for hypothesis testing.
#' @param alpha Significance level.
#' @param use_estm If `TRUE`, use [dcifer::ibdEstM()] instead of
#'   [dcifer::ibdPair()].
#' @param threads Number of parallel workers.
#' @param seed Random seed.
#' @param verbose Print parallel worker output.
#'
#' @return The relatedness tibble (also written to `btwn_host_rel_output`).
#' @export
dcifer_ibd_wrapper <- function(allele_table,
                               btwn_host_rel_output,
                               coi_table = NULL,
                               allele_freq_table = NULL,
                               specimen_name_col = "specimen_name",
                               target_name_col = "target_name",
                               target_value_col = "seq",
                               specimen_metadata = NULL,
                               pop_name_col = NULL,
                               rnull = 0,
                               alpha = 0.05,
                               use_estm = FALSE,
                               threads = 1L,
                               seed = 1L,
                               verbose = FALSE) {
  check_suggested_pkg("dcifer", "IBD relatedness via dcifer_ibd_wrapper()")
  check_suggested_pkgs(
    c("foreach", "iterators", "doParallel", "parallelly"),
    "parallel IBD estimation via dcifer_ibd_wrapper()"
  )

  if (is.null(allele_table) || is.null(btwn_host_rel_output)) {
    stop(
      "--allele_table and --btwn_host_rel_output are required",
      call. = FALSE
    )
  }
  set.seed(seed)

  allele_tbl <- create_dcifer_allele_table_input(
    allele_table,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col
  )
  dcifer_alleles <- dcifer::formatDat(
    allele_tbl,
    svar = "specimen_name",
    lvar = "target_name",
    avar = "target_value"
  )

  if (is.null(coi_table)) {
    coi <- dcifer::getCOI(dcifer_alleles)
    names(coi) <- stringr::str_split_i(names(coi), "\\.", 1)
  } else {
    coi <- create_dcifer_coi_input(
      coi_table,
      dcifer_alleles,
      specimen_name_col = specimen_name_col
    )
  }

  specimen_metadata_tbl <- NULL
  allele_freqs <- NULL
  if (is.null(allele_freq_table)) {
    if (is.null(pop_name_col)) {
      allele_freqs <- dcifer::calcAfreq(dcifer_alleles, coi, tol = 1e-5)
    } else {
      if (is.null(specimen_metadata)) {
        stop("pop_name_col provided but specimen_metadata missing", call. = FALSE)
      }
      specimen_metadata_tbl <- create_dcifer_specimen_metadata_input(
        specimen_metadata,
        coi,
        specimen_name_col = specimen_name_col,
        pop_name_col = pop_name_col
      )
    }
  } else {
    allele_freqs <- create_dcifer_allele_freq_input(
      allele_freq_table,
      dcifer_alleles,
      target_name_col = target_name_col,
      target_value_col = target_value_col
    )
  }

  if (is.null(pop_name_col)) {
    if (isTRUE(use_estm)) {
      dcifer_res <- run_ibdestm(
        dcifer_alleles,
        coi,
        allele_freqs,
        confint = TRUE,
        rnull = rnull,
        alpha = alpha,
        total_cores = threads,
        verbose = verbose
      )
    } else {
      dcifer_res <- run_ibdpair(
        dcifer_alleles,
        coi,
        allele_freqs,
        confint = TRUE,
        rnull = rnull,
        alpha = alpha,
        total_cores = threads,
        verbose = verbose
      )
    }
  } else {
    dcifer_res <- run_dcifer_bypop(
      allele_table,
      coi,
      specimen_metadata_tbl,
      specimen_name_col = specimen_name_col,
      target_name_col = target_name_col,
      target_value_col = target_value_col,
      use_estm = isTRUE(use_estm),
      confint = TRUE,
      rnull = rnull,
      alpha = alpha,
      total_cores = threads,
      verbose = verbose
    )
  }

  write_dcifer_ibd_output(dcifer_res, btwn_host_rel_output)
  dcifer_res
}
