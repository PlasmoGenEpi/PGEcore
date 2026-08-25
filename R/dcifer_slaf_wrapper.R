#' Read an allele table TSV for Dcifer wrappers
#'
#' @keywords internal
create_dcifer_allele_table_input <- function(allele_table_path,
                                             specimen_name_col = "specimen_name",
                                             target_name_col = "target_name",
                                             target_value_col = "seq") {
  allele_table <- readr::read_tsv(
    allele_table_path,
    col_types = do.call(
      readr::cols,
      c(
        stats::setNames(list(readr::col_character()), specimen_name_col),
        list(.default = readr::col_character())
      )
    ),
    progress = FALSE
  ) |>
    dplyr::select(dplyr::all_of(c(
      specimen_name_col, target_name_col, target_value_col
    ))) |>
    dplyr::rename(
      specimen_name = dplyr::all_of(specimen_name_col),
      target_name = dplyr::all_of(target_name_col),
      target_value = dplyr::all_of(target_value_col)
    )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(target_name),
    is.character(target_value),
    !is.na(specimen_name),
    !is.na(target_name),
    !is.na(target_value)
  )
  stop_on_validate_fails(allele_table, rules, "allele_table")
  allele_table
}

#' Read a COI table into a named vector aligned with a Dcifer allele list
#'
#' @keywords internal
create_dcifer_coi_input <- function(coi_path,
                                    allele_list,
                                    specimen_name_col = "specimen_name") {
  coi <- readr::read_tsv(
    coi_path,
    col_types = do.call(
      readr::cols,
      c(
        stats::setNames(list(readr::col_character()), specimen_name_col),
        list(.default = readr::col_character(), coi = readr::col_integer())
      )
    ),
    progress = FALSE
  ) |>
    dplyr::select(dplyr::all_of(specimen_name_col), "coi") |>
    dplyr::rename(specimen_name = dplyr::all_of(specimen_name_col))

  rules <- validate::validator(
    is.character(specimen_name),
    is.integer(coi),
    !is.na(specimen_name),
    !is.na(coi)
  )
  stop_on_validate_fails(coi, rules, "coi_table")

  specimen_inallele_notincoi <- setdiff(names(allele_list), coi$specimen_name)
  if (length(specimen_inallele_notincoi) > 0) {
    stop(
      "The following specimen IDs appear in the allele table and not in the ",
      "COI table: ",
      paste(specimen_inallele_notincoi, collapse = " "),
      call. = FALSE
    )
  }
  specimen_incoi_notinallele <- setdiff(coi$specimen_name, names(allele_list))
  if (length(specimen_incoi_notinallele) > 0) {
    warning(
      "The following specimen IDs appear in the allele table and not in the ",
      "COI table: ",
      paste(specimen_incoi_notinallele, collapse = " "),
      call. = FALSE
    )
  }

  coi <- tibble::tibble(specimen_name = names(allele_list)) |>
    dplyr::left_join(coi, by = "specimen_name")
  stats::setNames(coi$coi, coi$specimen_name)
}

#' Convert Dcifer allele-frequency lists to a tibble
#'
#' @keywords internal
prepare_dcifer_slaf_output <- function(allele_freqs_list,
                                       n_samp_per_target,
                                       target_name_col = "target_name",
                                       target_value_col = "seq") {
  onetargetaf_list2tib <- function(onetargetaf) {
    tibble::tibble(
      target_value = names(onetargetaf),
      freq = unname(onetargetaf)
    )
  }
  tibble::tibble(
    target_name = names(allele_freqs_list),
    alleles_freqs = unname(allele_freqs_list)
  ) |>
    dplyr::mutate(
      alleles_freqs = lapply(.data$alleles_freqs, onetargetaf_list2tib)
    ) |>
    tidyr::unnest("alleles_freqs") |>
    dplyr::left_join(n_samp_per_target, by = "target_name") |>
    dplyr::rename(
      !!target_name_col := "target_name",
      !!target_value_col := "target_value"
    )
}

#' Estimate single-locus allele frequencies with Dcifer
#'
#' File-oriented entry point used by the `dcifer_slaf_wrapper` CLI. The
#' **dcifer** package is an optional dependency (Suggests).
#'
#' @param allele_table Path to allele TSV.
#' @param slaf_output Path for SLAF TSV output.
#' @param coi_table Optional path to COI TSV.
#' @param specimen_name_col,target_name_col,target_value_col Column names.
#' @param tol,qstart Passed to `dcifer::calcAfreq()`.
#' @param coi_lrank Rank of the locus used by `dcifer::getCOI()` when
#'   `coi_table` is not supplied.
#'
#' @return The SLAF tibble (also written to `slaf_output`).
#' @export
dcifer_slaf_wrapper <- function(allele_table,
                                slaf_output,
                                coi_table = NULL,
                                specimen_name_col = "specimen_name",
                                target_name_col = "target_name",
                                target_value_col = "seq",
                                tol = 1e-04,
                                qstart = 0.5,
                                coi_lrank = 2L) {
  check_suggested_pkg(
    "dcifer",
    "single-locus allele frequencies via dcifer_slaf_wrapper()"
  )

  if (is.null(allele_table) || is.null(slaf_output)) {
    stop("--allele_table and --slaf_output are required", call. = FALSE)
  }
  if (!file.exists(allele_table)) {
    stop("allele_table file not found: ", allele_table, call. = FALSE)
  }

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
    coi <- dcifer::getCOI(dcifer_alleles, lrank = coi_lrank)
  } else {
    if (!identical(as.integer(coi_lrank), 2L)) {
      warning(
        "The --coi_lrank argument has been changed from the default, but this ",
        "will have no effect as --coi_table has also been specified",
        call. = FALSE
      )
    }
    coi <- create_dcifer_coi_input(
      coi_table,
      dcifer_alleles,
      specimen_name_col = specimen_name_col
    )
  }

  allele_freqs_list <- dcifer::calcAfreq(
    dcifer_alleles,
    coi,
    tol = tol,
    qstart = qstart
  )

  n_samp_per_target <- allele_tbl |>
    dplyr::group_by(.data$target_name) |>
    dplyr::summarize(
      sample_total = dplyr::n_distinct(.data$specimen_name),
      .groups = "drop"
    )

  slaf <- prepare_dcifer_slaf_output(
    allele_freqs_list,
    n_samp_per_target,
    target_name_col = target_name_col,
    target_value_col = target_value_col
  )
  readr::write_tsv(slaf, slaf_output)
  slaf
}
