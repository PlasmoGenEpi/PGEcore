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
  stop_on_validate_fails(loci_groups, rules, "loci_groups_input")

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
prepare_snpslice_af_output <- function(snpslice_res, loci_groups, use_mcmc) {
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
    use_map = !use_mcmc
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

  if (use_mcmc) {
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

#' Format SNP-Slice COI estimates
#'
#' @keywords internal
prepare_snpslice_coi_output <- function(snpslice_res,
                                        specimen_name_col,
                                        use_mcmc) {
  coi_tib <- snp.slicer::calculate_individual_coi(
    snpslice_res,
    use_map = !use_mcmc
  ) |>
    dplyr::select(-"host_index") |>
    dplyr::rename(!!specimen_name_col := "host_id", coi = "coi_estimate")
  if (!use_mcmc) {
    coi_tib <- dplyr::select(coi_tib, -"coi_sd", -"coi_lower", -"coi_upper")
  }
  coi_tib
}

#' Estimate multilocus allele frequency and COI with SNP-Slice
#'
#' File-oriented entry point used by the `snpslice_wrapper` CLI. Optional
#' **snp.slicer** and **variantstring** (1.x) must be installed separately.
#'
#' @param allele_table Path to allele TSV with counts.
#' @param loci_groups_input Path to loci-group TSV.
#' @param mlaf_output Path for multilocus allele-frequency TSV.
#' @param coi_output Path for COI TSV.
#' @param specimen_name_col,target_name_col,target_value_col,target_count_col
#'   Column names in `allele_table`.
#' @param loci_limit Optional cap on the number of loci.
#' @param model Observation model for SNP-Slice.
#' @param n_mcmc,burnin,alpha,rho,threshold,gap SNP-Slice MCMC settings.
#' @param use_mcmc_for_af_and_coi If `TRUE`, sample from MCMC instead of MAP.
#' @param verbose Verbose SNP-Slice output.
#' @param seed Random seed.
#'
#' @return Invisibly, a list with `mlaf` and `coi` tibbles.
#' @export
snpslice_wrapper <- function(allele_table,
                             loci_groups_input,
                             mlaf_output,
                             coi_output,
                             specimen_name_col = "specimen_name",
                             target_name_col = "aa_locus",
                             target_value_col = "aa",
                             target_count_col = "reads",
                             loci_limit = NULL,
                             model = "negative_binomial",
                             n_mcmc = 10000L,
                             burnin = NULL,
                             alpha = 2.6,
                             rho = 0.5,
                             threshold = 0.001,
                             gap = NULL,
                             use_mcmc_for_af_and_coi = FALSE,
                             verbose = FALSE,
                             seed = 1L) {
  check_suggested_pkg("snp.slicer", "SNP-Slice via snpslice_wrapper()")
  check_variantstring_v1("variant strings via snpslice_wrapper()")

  required <- list(
    allele_table = allele_table,
    loci_groups_input = loci_groups_input,
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

  set.seed(seed)
  allele_tbl <- create_snpslice_allele_table_input(
    allele_table,
    specimen_name_col = specimen_name_col,
    target_name_col = target_name_col,
    target_value_col = target_value_col,
    target_count_col = target_count_col
  )
  loci_groups <- create_snpslice_loci_group_input(
    loci_groups_input,
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
      loci_random <- sample(
        setdiff(allele_tbl$target_name, loci_oi),
        n_loci_select
      )
      loci_selected <- c(loci_oi, loci_random)
      allele_tbl <- dplyr::filter(allele_tbl, .data$target_name %in% loci_selected)
    }
  }

  snpslice_res <- snp.slicer::snp_slice(
    allele_tbl,
    model = model,
    n_mcmc = n_mcmc,
    burnin = burnin,
    alpha = alpha,
    rho = rho,
    threshold = threshold,
    gap = gap,
    store_mcmc = TRUE,
    verbose = verbose,
    specimen_id_col = "specimen_name",
    target_id_col = "target_name",
    target_value_col = "target_value",
    target_count_col = "target_count"
  )

  mlaf <- prepare_snpslice_af_output(
    snpslice_res,
    loci_groups,
    use_mcmc_for_af_and_coi
  )
  readr::write_tsv(mlaf, mlaf_output)
  coi <- prepare_snpslice_coi_output(
    snpslice_res,
    specimen_name_col,
    use_mcmc_for_af_and_coi
  )
  readr::write_tsv(coi, coi_output)
  invisible(list(mlaf = mlaf, coi = coi))
}
