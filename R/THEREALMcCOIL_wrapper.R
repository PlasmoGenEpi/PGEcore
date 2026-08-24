#' Read SNP call data for THEREALMcCOIL
#'
#' @param snp_calls_input Path to an independent, collapsed SNP call TSV.
#' @return A data frame with `specimen_name`, `snp_name`, `seq_base`, `reads`.
#' @keywords internal
read_and_preprocess_snp_call <- function(snp_calls_input) {
  required_cols <- c("specimen_name", "snp_name", "reads", "seq_base")
  df_snp_call <- readr::read_tsv(
    snp_calls_input,
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

#' Format SNP calls for the McCOIL proportional model
#'
#' @param df Output of `read_and_preprocess_snp_call()`.
#' @return A list with `a1` and `a2` read-count matrices.
#' @keywords internal
prep_input_prop <- function(df) {
  allele_map <- df |>
    dplyr::arrange(.data$snp_name, .data$seq_base) |>
    dplyr::distinct(.data$snp_name, .data$seq_base) |>
    dplyr::mutate(allele_idx = rep_len(c(1, 2), length.out = dplyr::n()))
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

  list(a1 = df_allele1, a2 = df_allele2)
}

#' Run McCOIL categorical or proportional MCMC into `work_dir`
#'
#' @param df Preprocessed SNP calls.
#' @param model `"categorical"` or `"proportional"`.
#' @param work_dir Directory for McCOIL temp traces (typically `tempdir()`).
#' @param output Base filename written under `work_dir`.
#' @keywords internal
call_mccoil <- function(df,
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
                        work_dir,
                        output = "McCOIL_out.txt") {
  if (!model %in% c("categorical", "proportional")) {
    stop("--model must be one of categorical|proportional", call. = FALSE)
  }
  if (!err_method %in% c(1, 3)) {
    stop("--err_method must be one of 1|3", call. = FALSE)
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  if (model == "categorical") {
    mccoil_cat_input <- prep_input_categorical(df)
    run_mccoil_categorical(
      mccoil_cat_input,
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
    mccoil_prop_input <- prep_input_prop(df)
    run_mccoil_proportional(
      mccoil_prop_input$a1,
      mccoil_prop_input$a2,
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
  invisible(file.path(work_dir, paste0(output, "_summary.txt")))
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
#' File-oriented entry point used by the `THEREALMcCOIL_wrapper` CLI. Compiled
#' C MCMC routines (`McCOIL_categorical`, `McCOIL_prop`) are linked at package
#' install time.
#'
#' @param snp_calls_input TSV of SNP calls with at least `specimen_name`,
#'   `snp_name`, `pos`, `seq_base`, `reads` (and typically `target_name`, `he`).
#' @param slaf_output Output TSV of allele frequencies (`variant`, `freq`).
#' @param coi_output Output TSV of COI estimates (`specimen_name`, `coi`).
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
#'
#' @return A list with `slaf` and `coi` (invisibly after writing outputs).
#' @export
THEREALMcCOIL_wrapper <- function(snp_calls_input,
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
                                  err_method = 1L) {
  if (mccoil_blank(snp_calls_input)) {
    stop("--snp_calls_input must be set", call. = FALSE)
  }
  if (mccoil_blank(slaf_output)) {
    stop("--slaf_output must be set", call. = FALSE)
  }
  if (mccoil_blank(coi_output)) {
    stop("--coi_output must be set", call. = FALSE)
  }

  df <- read_and_preprocess_snp_call(snp_calls_input)
  work_dir <- tempfile("McCOIL_")
  dir.create(work_dir)
  on.exit(clean_up_mccoil(work_dir), add = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  summary_path <- call_mccoil(
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
    work_dir = work_dir,
    output = "McCOIL_out.txt"
  )
  df_formated <- format_mccoil_output(summary_path)
  write_mccoil_output(df_formated, slaf_output, coi_output)
  invisible(df_formated)
}
