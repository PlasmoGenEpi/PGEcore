#' Parse a comma-separated name list or a one-column TSV of names
#'
#' @param value Character scalar, or `NULL`.
#' @return Character vector of names (possibly empty).
#' @keywords internal
parse_name_list_arg <- function(value) {
  if (is.null(value)) {
    return(character(0))
  }
  if (length(value) > 1L) {
    return(as.character(value))
  }
  if (!nzchar(value)) {
    return(character(0))
  }
  if (file.exists(value)) {
    as.character(
      readr::read_tsv(value, col_names = FALSE, show_col_types = FALSE)[[1]]
    )
  } else {
    as.character(unlist(strsplit(value, split = ",", fixed = TRUE)))
  }
}

#' Filter an SNP table by optional target and specimen selections
#'
#' @param snp_data SNP call tibble.
#' @param select_target_names Character vector of target names (empty = no filter).
#' @param select_specimen_names Character vector of specimen names (empty = no filter).
#' @return Filtered `snp_data`.
#' @keywords internal
filter_snp_table_for_optional_subselecting <- function(snp_data,
                                                      select_target_names = character(0),
                                                      select_specimen_names = character(0)) {
  if (length(select_specimen_names) > 0) {
    snp_data <- dplyr::filter(
      snp_data,
      .data$specimen_name %in% select_specimen_names
    )
  }
  if (length(select_target_names) > 0) {
    snp_data <- dplyr::filter(
      snp_data,
      .data$target_name %in% select_target_names
    )
  }
  snp_data
}

#' Stop when requested specimen/target names are missing from SNP data
#'
#' @param snp_data SNP call tibble.
#' @param select_target_names Character vector of requested targets.
#' @param select_specimen_names Character vector of requested specimens.
#' @param snp_table_fnp Path label used in messages.
#' @return Invisibly `TRUE` if all requested names are present.
#' @keywords internal
stop_on_missing_subselections <- function(snp_data,
                                          select_target_names,
                                          select_specimen_names,
                                          snp_table_fnp) {
  warns <- character(0)
  if (length(select_specimen_names) > 0) {
    missing_sel_specs <- setdiff(
      select_specimen_names,
      unique(snp_data$specimen_name)
    )
    if (length(missing_sel_specs) > 0) {
      warns <- c(
        warns,
        paste0(
          "supplied --select_specimen_names but the following specimen_names ",
          "are missing from ", snp_table_fnp, "\n",
          paste(missing_sel_specs, collapse = ",")
        )
      )
    }
  }
  if (length(select_target_names) > 0) {
    missing_sel_tars <- setdiff(
      select_target_names,
      unique(snp_data$target_name)
    )
    if (length(missing_sel_tars) > 0) {
      warns <- c(
        warns,
        paste0(
          "supplied --select_target_names but the following target_names ",
          "are missing from ", snp_table_fnp, "\n",
          paste(missing_sel_tars, collapse = ",")
        )
      )
    }
  }
  if (length(warns) > 0) {
    stop(paste0("\n", paste(warns, collapse = "\n")), call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate SNP table column types for independent-SNP filtering
#'
#' @param snp_table SNP call tibble.
#' @return Invisibly `TRUE` if valid.
#' @keywords internal
validate_snp_table_column_types <- function(snp_table) {
  snp_table_rules <- validate::validator(
    is.character(specimen_name),
    is.numeric(reads),
    is.character(target_name),
    is.character(ref_base),
    is.character(seq_base),
    is.character(chrom),
    is.character(snp_name),
    is.numeric(pos),
    is.logical(is_biallelic),
    !is.na(specimen_name),
    !is.na(reads),
    !is.na(target_name),
    !is.na(ref_base),
    !is.na(seq_base),
    !is.na(chrom),
    !is.na(snp_name),
    !is.na(pos),
    !is.na(is_biallelic)
  )
  stop_on_validate_fails(snp_table, snp_table_rules, "snp_table")
}

#' Core greedy filter on expected heterozygosity and distance
#'
#' @keywords internal
filter_highest_diversity_snps_core <- function(snp_table_in,
                                               mindist_between_snps = 10000,
                                               only_biallelic = FALSE,
                                               only_informative = FALSE) {
  if (only_biallelic) {
    snp_table_in <- dplyr::filter(snp_table_in, .data$is_biallelic)
  }

  snp_table_counts_per_specimen_name <- snp_table_in |>
    dplyr::group_by(
      .data$specimen_name, .data$chrom, .data$pos, .data$snp_name,
      .data$ref_base, .data$seq_base
    ) |>
    dplyr::mutate(n = dplyr::n()) |>
    dplyr::ungroup()

  snp_table_counts_per_specimen_name_multi <- dplyr::filter(
    snp_table_counts_per_specimen_name,
    .data$n > 1
  )
  if (nrow(snp_table_counts_per_specimen_name_multi) > 0) {
    stop(
      "the following snps were found to have multiple calls per specimen_name ",
      "for the same seq_base, make sure calls are collapsed:\n",
      paste(
        unique(snp_table_counts_per_specimen_name_multi$snp_name),
        collapse = ","
      ),
      call. = FALSE
    )
  }

  snp_table_freqs <- snp_table_in |>
    dplyr::group_by(
      .data$chrom, .data$pos, .data$snp_name, .data$ref_base, .data$seq_base
    ) |>
    dplyr::summarise(allele_count = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(.data$chrom, .data$pos, .data$snp_name, .data$ref_base) |>
    dplyr::mutate(
      allele_total = sum(.data$allele_count),
      allele_freq = .data$allele_count / .data$allele_total
    )

  snp_table_he <- snp_table_freqs |>
    dplyr::group_by(.data$chrom, .data$pos, .data$snp_name, .data$ref_base) |>
    dplyr::summarise(he = 1 - sum(.data$allele_freq^2), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$he))

  snp_table_he$keep <- TRUE
  if (nrow(snp_table_he) >= 2) {
    for (row in seq.int(2L, nrow(snp_table_he))) {
      for (filt_row in seq_len(row - 1L)) {
        if (
          snp_table_he$keep[filt_row] &&
            snp_table_he$chrom[filt_row] == snp_table_he$chrom[row] &&
            abs(snp_table_he$pos[filt_row] - snp_table_he$pos[row]) <
              mindist_between_snps
        ) {
          snp_table_he$keep[row] <- FALSE
          break
        }
      }
    }
  }

  snp_table_he_filt <- snp_table_he |>
    dplyr::filter(.data$keep) |>
    dplyr::mutate(snp_id = paste0(.data$chrom, "-", .data$pos))

  if (only_informative) {
    snp_table_he_filt <- dplyr::filter(snp_table_he_filt, .data$he > 0)
  }

  snp_table_in <- snp_table_in |>
    dplyr::left_join(
      dplyr::select(snp_table_he, -"keep"),
      by = c("chrom", "pos", "snp_name", "ref_base")
    )

  snp_table_in |>
    dplyr::filter(
      paste0(.data$chrom, "-", .data$pos) %in% snp_table_he_filt$snp_id
    )
}

#' Filter SNPs to highest-diversity loci spaced by a minimum distance
#'
#' Ranks SNPs by expected heterozygosity, then greedily keeps loci that are at
#' least `mindist_between_snps` apart on the same chromosome.
#'
#' @param snp_table_in Path to a SNP TSV, or a data frame, with columns
#'   `specimen_name`, `target_name`, `chrom`, `pos`, `snp_name`, `ref_base`,
#'   `seq_base`, `reads`, and `is_biallelic`.
#' @param snp_table_out Optional output TSV path.
#' @param mindist_between_snps Minimum distance between kept SNPs (default
#'   `10000`).
#' @param select_target_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of targets to keep.
#' @param select_specimen_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of specimens to keep.
#' @param overwrite If `FALSE` (default), refuse to overwrite `snp_table_out`.
#' @param only_biallelic If `TRUE`, restrict to rows where `is_biallelic` is
#'   `TRUE`.
#' @param only_informative If `TRUE`, drop SNPs with expected heterozygosity
#'   of zero.
#'
#' @return Filtered SNP table including an `he` column.
#'
#' @export
filter_to_highest_diversity_independent_snp_call <- function(snp_table_in,
                                                             snp_table_out = NULL,
                                                             mindist_between_snps = 10000,
                                                             select_target_names = NULL,
                                                             select_specimen_names = NULL,
                                                             overwrite = FALSE,
                                                             only_biallelic = FALSE,
                                                             only_informative = FALSE) {
  options(dplyr.summarise.inform = FALSE)

  input_label <- "snp_table_in"
  if (is.character(snp_table_in) && length(snp_table_in) == 1L) {
    if (!file.exists(snp_table_in)) {
      stop(snp_table_in, " does not exist", call. = FALSE)
    }
    input_label <- snp_table_in
    snp_tbl <- readr::read_tsv(
      snp_table_in,
      col_types = readr::cols(specimen_name = readr::col_character()),
      show_col_types = FALSE
    )
  } else if (is.data.frame(snp_table_in)) {
    snp_tbl <- tibble::as_tibble(snp_table_in)
  } else {
    stop("`snp_table_in` must be a file path or a data frame.", call. = FALSE)
  }

  validate_required_columns(
    snp_tbl,
    c(
      "specimen_name", "target_name", "chrom", "pos", "snp_name",
      "ref_base", "seq_base", "reads", "is_biallelic"
    ),
    input_label
  )
  validate_snp_table_column_types(snp_tbl)

  select_targets <- parse_name_list_arg(select_target_names)
  select_specs <- parse_name_list_arg(select_specimen_names)
  stop_on_missing_subselections(
    snp_tbl,
    select_targets,
    select_specs,
    input_label
  )
  snp_tbl <- filter_snp_table_for_optional_subselecting(
    snp_tbl,
    select_targets,
    select_specs
  )

  stop_if_output_exists(snp_table_out, overwrite)

  out <- filter_highest_diversity_snps_core(
    snp_tbl,
    mindist_between_snps = mindist_between_snps,
    only_biallelic = only_biallelic,
    only_informative = only_informative
  )

  if (!is.null(snp_table_out)) {
    readr::write_tsv(out, snp_table_out)
  }
  out
}
