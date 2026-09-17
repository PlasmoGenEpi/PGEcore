#' Treat NULL or blank CLI strings as missing
#'
#' @noRd
idm_blank <- function(x) {
  is.null(x) || (is.character(x) && !nzchar(x))
}

#' Prepare allele-table input for the Incomplete Data Model
#'
#' @param allele_table Path to allele TSV.
#' @return Tibble with `specimen_name`, `locus`, `variants`.
#' @keywords internal
prepare_input_4_allele_table <- function(allele_table) {
  df <- readr::read_tsv(
    allele_table,
    col_types = readr::cols(
      .default = readr::col_character(),
      specimen_name = readr::col_character(),
      reads = readr::col_integer()
    )
  )
  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(seq),
    is.integer(reads),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(seq),
    !is.na(reads)
  )
  stop_on_validate_fails(df, rules, "allele_table")
  df |>
    dplyr::mutate(
      locus = .data$target_name,
      variants = stringr::str_c(.data$target_name, .data$seq, sep = ":")
    ) |>
    dplyr::select("specimen_name", "locus", "variants") |>
    dplyr::distinct(.data$specimen_name, .data$locus, .data$variants)
}

#' Prepare amino-acid-call input for the Incomplete Data Model
#'
#' @param aa_calls Path to amino-acid TSV.
#' @return Tibble with `specimen_name`, `locus`, `variants`.
#' @keywords internal
prepare_input_4_aa_calls <- function(aa_calls) {
  df <- readr::read_tsv(
    aa_calls,
    col_types = readr::cols(
      .default = readr::col_character(),
      specimen_name = readr::col_character(),
      aa_position = readr::col_integer()
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
  stop_on_validate_fails(df, rules, "aa_calls")
  df |>
    dplyr::mutate(
      locus = stringr::str_c(.data$gene_id, .data$aa_position, sep = ":"),
      variants = stringr::str_c(.data$gene_id, .data$aa_position, .data$aa, sep = ":")
    ) |>
    dplyr::select("specimen_name", "locus", "variants")
}

#' Solve one locus with the vendored MLE
#'
#' @param tbl Two-column table (`specimen_name`, `variants`) for a single locus.
#' @param model `"IDM"` or `"OM"` — the value handed to the vendored `MLE()`.
#' @param continuity Continuity correction applied to lineage prevalence; `0`
#'   is the published estimator.
#' @param lambda_initial Initial lambda for the numerical solver.
#' @param eps_initial Initial epsilon for the numerical solver.
#' @return List with `variants` and `freq`, one entry per variant at the locus.
#' @keywords internal
idm_locus_mle <- function(tbl, model, lambda_initial, eps_initial,
                          continuity = 0) {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  utils::write.table(tbl, tmp, sep = "\t", row.names = FALSE)
  dat <- .idm_vendor$DatImp(tmp)
  nk <- .idm_vendor$Nk(dat)
  mle_res <- .idm_vendor$MLE(
    nk[[1]],
    nk[[2]],
    nk[[3]],
    model = model,
    lambda_initial = lambda_initial,
    eps_initial = eps_initial,
    continuity = continuity
  )
  n_uniq_variants <- length(nk$N_k)
  locus_freq <- mle_res$`lineage frequencies`
  if (length(locus_freq) < n_uniq_variants) {
    locus_freq <- rep(NA, n_uniq_variants)
  }
  list(variants = colnames(nk$N_k), freq = locus_freq)
}

#' Run IDM/OM MLE independently at each locus
#'
#' @param df Formatted input from `prepare_input_4_allele_table()` or
#'   `prepare_input_4_aa_calls()`.
#' @param model `"IDM"`, `"OM"`, `"IDM_OM"`, `"OM_CC"`, or `"IDM_OM_CC"`.
#' @param lambda_initial Initial lambda for the numerical solver.
#' @param eps_initial Initial epsilon for the numerical solver.
#' @return Tibble with `variant` and `freq`.
#' @keywords internal
run_idm_mle_across_loci <- function(df,
                                    model = "IDM",
                                    lambda_initial = 1.0,
                                    eps_initial = 0.1) {
  variants_array <- c()
  freq_array <- c()

  # An "IDM_OM" model runs the IDM and falls back to the OM at any locus the IDM
  # leaves unsolved. A "_CC" model applies a continuity correction to lineage
  # prevalence, keeping a lineage seen in every sample off the boundary where the
  # MLE does not exist; it applies to whichever OM solve the model reaches -- the
  # only solve for "OM_CC", the fallback solve for "IDM_OM_CC" -- while the IDM
  # solve itself always runs uncorrected. The vendored MLE() only knows "IDM"
  # and "OM".
  fallback_to_om <- model %in% c("IDM_OM", "IDM_OM_CC")
  om_continuity <- if (grepl("_CC$", model)) 0.5 else 0
  primary_continuity <- if (fallback_to_om) 0 else om_continuity
  locus_model <- if (fallback_to_om) "IDM" else sub("_CC$", "", model)

  # n_0, the count of specimens with no lineage detected, is what separates the
  # IDM from the OM: with n_0 == 0 the IDM branch of MLE1() falls through to
  # MLE_OM(). Nk() derives n_0 as N - Nplus from the table it is handed, so a
  # specimen absent at a locus has to appear there carrying NA or the
  # incomplete-data model never runs.
  #
  # Only the IDM accepts that. MLE() rejects `sum(N_k) < N` outright for the OM,
  # whose contract is complete data, so the OM keeps seeing observed calls only.
  specimens <- if (identical(locus_model, "IDM")) {
    dplyr::distinct(df, .data$specimen_name)
  } else {
    NULL
  }

  for (l in dplyr::distinct(df, .data$locus)$locus) {
    observed <- df |>
      dplyr::filter(.data$locus == l) |>
      dplyr::select("specimen_name", "variants")
    tmp_df <- if (is.null(specimens)) {
      observed
    } else {
      dplyr::left_join(
        specimens, observed,
        by = "specimen_name",
        relationship = "many-to-many"
      )
    }
    fit <- idm_locus_mle(
      tmp_df, locus_model, lambda_initial, eps_initial, primary_continuity
    )
    # A locus the IDM cannot solve comes back all-NA. Replacing the whole locus
    # rather than the individual NA entries keeps the returned frequencies a
    # single coherent solution for that locus.
    if (fallback_to_om && anyNA(fit$freq)) {
      fit <- idm_locus_mle(
        observed, "OM", lambda_initial, eps_initial, om_continuity
      )
    }
    freq_array <- c(freq_array, fit$freq)
    variants_array <- c(variants_array, fit$variants)
  }
  tibble::tibble(variant = variants_array, freq = freq_array)
}

#' Write IDM SLAF output
#'
#' @param res Result table from `run_idm_mle_across_loci()`.
#' @param slaf_output Output TSV path.
#' @param allele_table If `TRUE`, split `variant` into `target_name` and `seq`.
#' @keywords internal
write_idm_output <- function(res, slaf_output, allele_table = FALSE) {
  if (allele_table) {
    res |>
      tidyr::separate_wider_delim(
        "variant",
        ":",
        names = c("target_name", "seq")
      ) |>
      readr::write_tsv(slaf_output)
  } else {
    readr::write_tsv(res, slaf_output)
  }
}

#' Estimate single-locus allele frequencies with the Incomplete Data Model
#'
#' Estimates single-locus allele frequencies with the Incomplete Data Model
#' (or original model). Provide exactly one of `allele_table` or `aa_calls`.
#' Vendored MLE code is from Hashemi & Schneider (2024). Requires **Rmpfr** and
#' **openxlsx** (Suggests).
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele table TSV, or `""` / `NULL` if using
#'   `aa_calls`. See `vignette("input-formats", package = "PGEcore")`.
#' - **`aa_calls`**: Amino-acid calls TSV, or `""` / `NULL` if using
#'   `allele_table`.
#'
#' ## Outputs
#'
#' - **`slaf_output`**: Single-locus allele frequencies. Allele-table input is
#'   written as `target_name`, `seq`, `freq`; AA-call input as `variant`,
#'   `freq`.
#'
#' ## Running
#'
#' ```r
#' IDM_wrapper(
#'   allele_table = "allele_table.tsv",
#'   slaf_output = "slaf.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/IDM_wrapper \
#'   --allele_table allele_table.tsv \
#'   --slaf_output slaf.tsv
#' ```
#'
#' Requires **Rmpfr** and **openxlsx** (Suggests).
#'
#' @param allele_table Path to allele table TSV, or `""` / `NULL` if using
#'   amino-acid calls. See *Inputs*.
#' @param aa_calls Path to amino-acid calls TSV, or `""` / `NULL` if using an
#'   allele table. See *Inputs*.
#' @param slaf_output Output TSV path. See *Outputs*.
#' @param model `"IDM"` (incomplete-data model), `"OM"` (original model), or
#'   `"IDM_OM"` (IDM, falling back to the OM at any locus the IDM leaves
#'   unsolved), `"OM_CC"` (OM with a continuity correction on lineage
#'   prevalence), or `"IDM_OM_CC"` (IDM, falling back to `"OM_CC"`).
#' @param lambda_initial Initial lambda for the numerical iteration.
#' @param eps_initial Initial epsilon for the numerical iteration.
#'
#' @return The SLAF tibble (invisibly after writing `slaf_output`).
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
IDM_wrapper <- function(allele_table = "",
                        aa_calls = "",
                        slaf_output,
                        model = "IDM",
                        lambda_initial = 1.0,
                        eps_initial = 0.1) {
  if (idm_blank(slaf_output)) {
    stop("Missing required arguments: --slaf_output", call. = FALSE)
  }
  n_inputs <- as.integer(!idm_blank(allele_table)) +
    as.integer(!idm_blank(aa_calls))
  if (n_inputs != 1L) {
    stop(
      "One and only one of the args --allele_table and ",
      "--aa_calls must be provided.",
      call. = FALSE
    )
  }
  if (!model %in% c("IDM", "OM", "IDM_OM", "OM_CC", "IDM_OM_CC")) {
    stop(
      "--model must be one of IDM | OM | IDM_OM | OM_CC | IDM_OM_CC",
      call. = FALSE
    )
  }
  check_suggested_pkg("Rmpfr", "Incomplete Data Model MLE via IDM_wrapper()")
  check_suggested_pkg("openxlsx", "Incomplete Data Model data import via IDM_wrapper()")

  if (!idm_blank(aa_calls)) {
    df <- prepare_input_4_aa_calls(aa_calls)
  } else {
    df <- prepare_input_4_allele_table(allele_table)
  }
  res <- run_idm_mle_across_loci(
    df,
    model,
    lambda_initial,
    eps_initial
  )
  write_idm_output(
    res,
    slaf_output,
    allele_table = idm_blank(aa_calls)
  )
  invisible(res)
}
