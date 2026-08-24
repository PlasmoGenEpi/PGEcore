#' Read amino acid calls for naive multilocus prev/freq
#'
#' @param aa_table Path to a TSV of amino acid calls.
#' @return A tibble including `n_aa` (allele count per specimen/locus).
#' @keywords internal
read_mlp_naive_aa_table <- function(aa_table) {
  stopifnot(is.character(aa_table), length(aa_table) == 1L)
  if (!file.exists(aa_table)) {
    stop(aa_table, " does not exist", call. = FALSE)
  }

  # Match legacy read.table() integer types used by validate rules.
  df_aa <- utils::read.table(
    aa_table,
    header = TRUE,
    colClasses = c(specimen_name = "character")
  )

  validate_required_columns(
    df_aa,
    c("specimen_name", "gene", "gene_id", "aa_position", "reads", "aa"),
    "aa_table"
  )

  rules <- validate::validator(
    is.character(specimen_name),
    is.character(gene),
    is.character(gene_id),
    is.integer(aa_position),
    is.integer(reads),
    is.character(aa),
    !is.na(specimen_name),
    !is.na(gene),
    !is.na(gene_id),
    !is.na(aa_position),
    !is.na(reads),
    !is.na(aa)
  )
  stop_on_validate_fails(df_aa, rules, "aa_table")

  df_aa |>
    dplyr::group_by(.data$specimen_name, .data$gene_id, .data$aa_position) |>
    dplyr::mutate(n_aa = dplyr::n()) |>
    dplyr::ungroup()
}

#' Read loci-group definitions for naive multilocus prev/freq
#'
#' @param loci_groups_path Path to a TSV with `group_id`, `gene_id`, `aa_position`.
#' @return A tibble of loci groups.
#' @keywords internal
read_mlp_naive_loci_groups <- function(loci_groups_path) {
  stopifnot(is.character(loci_groups_path), length(loci_groups_path) == 1L)
  if (!file.exists(loci_groups_path)) {
    stop(loci_groups_path, " does not exist", call. = FALSE)
  }

  loci_groups <- readr::read_tsv(
    loci_groups_path,
    col_types = readr::cols(
      .default = readr::col_character(),
      aa_position = readr::col_integer()
    ),
    progress = FALSE
  )

  validate_required_columns(
    loci_groups,
    c("group_id", "gene_id", "aa_position"),
    "loci_groups_input"
  )

  rules <- validate::validator(
    is.character(group_id),
    is.character(gene_id),
    is.integer(aa_position),
    !is.na(group_id),
    !is.na(gene_id),
    !is.na(aa_position)
  )
  stop_on_validate_fails(loci_groups, rules, "loci_groups_input")
  loci_groups
}

#' Recalculate single-locus prev/freq from multilocus calls (WSAF-weighted)
#'
#' @param multilocus_calls Table with `group_id`, `specimen_name`, `variant`, `wsaf`.
#' @return Summarized single-locus prev/freq per group.
#' @keywords internal
generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop <- function(multilocus_calls) {
  rules <- validate::validator(
    is.character(group_id),
    is.character(specimen_name),
    is.character(variant),
    is.numeric(wsaf),
    !is.na(group_id),
    !is.na(specimen_name),
    !is.na(variant),
    !is.na(wsaf)
  )
  stop_on_validate_fails(
    multilocus_calls,
    rules,
    "multilocus_calls for generate_single_locus_prev_freq_from_multilocus_groups"
  )

  multilocus_calls_mod <- multilocus_calls |>
    dplyr::group_by(.data$group_id) |>
    dplyr::mutate(variant_split = strsplit(.data$variant, split = ";")) |>
    tidyr::unnest("variant_split") |>
    tidyr::separate(
      "variant_split",
      into = c("gene_id", "aa_position", "aa"),
      sep = ":"
    ) |>
    dplyr::mutate(
      aa_position = strsplit(.data$aa_position, split = "_"),
      aa = strsplit(.data$aa, split = "_")
    ) |>
    tidyr::unnest(c("aa_position", "aa")) |>
    dplyr::mutate(aa_position = as.numeric(.data$aa_position)) |>
    dplyr::arrange(.data$group_id, .data$gene_id, .data$aa_position, .data$aa)

  multilocus_calls_mod |>
    dplyr::group_by(.data$group_id, .data$gene_id, .data$aa_position) |>
    dplyr::mutate(
      sample_total = dplyr::n_distinct(.data$specimen_name),
      wsaf_total = sum(.data$wsaf)
    ) |>
    dplyr::group_by(
      .data$group_id, .data$gene_id, .data$aa_position, .data$aa,
      .data$sample_total, .data$wsaf_total
    ) |>
    dplyr::summarise(
      wsaf_sum = sum(.data$wsaf),
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(
      freq = .data$wsaf_sum / .data$wsaf_total,
      prev = .data$sample_count / .data$sample_total
    ) |>
    tidyr::unite("variant", "gene_id", "aa_position", "aa", sep = ":") |>
    dplyr::select(-"wsaf_sum", -"wsaf_total")
}

#' Recalculate single-locus prev/freq from multilocus calls (presence/absence)
#'
#' @param multilocus_calls Table with `group_id`, `specimen_name`, `variant`.
#' @return Summarized single-locus prev/freq per group.
#' @keywords internal
generate_single_locus_prev_freq_from_multilocus_groups_presence_absence <- function(multilocus_calls) {
  rules <- validate::validator(
    is.character(group_id),
    is.character(specimen_name),
    is.character(variant),
    !is.na(group_id),
    !is.na(specimen_name),
    !is.na(variant)
  )
  stop_on_validate_fails(
    multilocus_calls,
    rules,
    "multilocus_calls for generate_single_locus_prev_freq_from_multilocus_groups"
  )

  multilocus_calls_mod <- multilocus_calls |>
    dplyr::group_by(.data$group_id) |>
    dplyr::mutate(variant_split = strsplit(.data$variant, split = ";")) |>
    tidyr::unnest("variant_split") |>
    tidyr::separate(
      "variant_split",
      into = c("gene_id", "aa_position", "aa"),
      sep = ":"
    ) |>
    dplyr::mutate(
      aa_position = strsplit(.data$aa_position, split = "_"),
      aa = strsplit(.data$aa, split = "_")
    ) |>
    tidyr::unnest(c("aa_position", "aa")) |>
    dplyr::mutate(aa_position = as.numeric(.data$aa_position)) |>
    dplyr::arrange(.data$group_id, .data$gene_id, .data$aa_position, .data$aa)

  multilocus_calls_mod |>
    dplyr::group_by(.data$group_id, .data$gene_id, .data$aa_position) |>
    dplyr::mutate(
      sample_total = dplyr::n_distinct(.data$specimen_name),
      allele_total = dplyr::n()
    ) |>
    dplyr::group_by(
      .data$group_id, .data$gene_id, .data$aa_position, .data$aa,
      .data$sample_total, .data$allele_total
    ) |>
    dplyr::summarise(
      allele_count = dplyr::n(),
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(
      freq = .data$allele_count / .data$allele_total,
      prev = .data$sample_count / .data$sample_total
    ) |>
    tidyr::unite("variant", "gene_id", "aa_position", "aa", sep = ":")
}

#' Calculate multilocus prev/freq weighted by WSAF
#'
#' @param multilocus_calls Table with `specimen_name`, `variant`, `wsaf`.
#' @return A tibble with `prev` and `freq`.
#' @keywords internal
calculate_multilocus_af_prev_wsaf_prop <- function(multilocus_calls) {
  rules <- validate::validator(
    is.character(specimen_name),
    is.character(variant),
    is.numeric(wsaf),
    !is.na(specimen_name),
    !is.na(variant),
    !is.na(wsaf)
  )
  stop_on_validate_fails(
    multilocus_calls,
    rules,
    "multilocus_calls for calculate_multilocus_af_prev_presence_absence"
  )

  multilocus_calls |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sample_total = dplyr::n_distinct(.data$specimen_name),
      wsaf_total = sum(.data$wsaf)
    ) |>
    dplyr::group_by(.data$variant, .data$sample_total, .data$wsaf_total) |>
    dplyr::summarise(
      wsaf_sum = sum(.data$wsaf),
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(
      prev = .data$sample_count / .data$sample_total,
      freq = .data$wsaf_sum / .data$wsaf_total
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-"wsaf_sum", -"wsaf_total")
}

#' Calculate multilocus prev/freq by presence/absence
#'
#' @param multilocus_calls Table with `specimen_name`, `variant`.
#' @return A tibble with `prev` and `freq`.
#' @keywords internal
calculate_multilocus_af_prev_presence_absence <- function(multilocus_calls) {
  rules <- validate::validator(
    is.character(specimen_name),
    is.character(variant),
    !is.na(specimen_name),
    !is.na(variant)
  )
  stop_on_validate_fails(
    multilocus_calls,
    rules,
    "multilocus_calls for calculate_multilocus_af_prev_presence_absence"
  )

  multilocus_calls |>
    dplyr::ungroup() |>
    dplyr::mutate(
      sample_total = dplyr::n_distinct(.data$specimen_name),
      allele_total = dplyr::n()
    ) |>
    dplyr::group_by(.data$variant, .data$sample_total, .data$allele_total) |>
    dplyr::summarise(
      allele_count = dplyr::n(),
      sample_count = dplyr::n_distinct(.data$specimen_name)
    ) |>
    dplyr::mutate(
      prev = .data$sample_count / .data$sample_total,
      freq = .data$allele_count / .data$allele_total
    ) |>
    dplyr::ungroup()
}

#' Build naive multilocus haplotype calls for one loci group
#'
#' @param aa_table Amino acid calls with `n_aa`.
#' @param group_df Loci-group rows including `loci_in_group`.
#' @param wsaf_cut_off Dominant-allele WSAF threshold.
#' @return A tibble of specimen haplotypes (`specimen_name`, `variant`, `wsaf`).
#' @keywords internal
build_naive_multilocus_calls_for_group <- function(aa_table, group_df, wsaf_cut_off) {
  aa_table_group <- aa_table |>
    dplyr::inner_join(group_df, by = c("gene_id", "aa_position")) |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::mutate(
      loci_called = dplyr::n_distinct(paste0(.data$gene_id, "-", .data$aa_position))
    ) |>
    dplyr::filter(.data$loci_called == .data$loci_in_group) |>
    dplyr::ungroup()

  aa_table_group_only_1_variable <- aa_table_group |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::filter(sum(.data$n_aa == 1) == (.data$loci_in_group - 1))

  if (nrow(aa_table_group_only_1_variable) > 0) {
    aa_table_group_only_1_variable_filt_variable <- aa_table_group_only_1_variable |>
      dplyr::filter(.data$n_aa != 1) |>
      dplyr::group_by(.data$specimen_name, .data$gene, .data$gene_id, .data$aa_position) |>
      dplyr::mutate(total_reads = sum(.data$reads)) |>
      dplyr::mutate(wsaf = .data$reads / .data$total_reads) |>
      dplyr::group_by(.data$specimen_name) |>
      dplyr::mutate(within_sample_hap = dplyr::row_number())

    aa_table_group_only_1_variable_filt_invariable <- aa_table_group_only_1_variable |>
      dplyr::filter(.data$n_aa == 1) |>
      dplyr::ungroup() |>
      dplyr::select(
        "specimen_name", "gene_id", "aa_position", "aa", "group_id"
      ) |>
      dplyr::left_join(
        aa_table_group_only_1_variable_filt_variable |>
          dplyr::ungroup() |>
          dplyr::group_by(.data$specimen_name) |>
          dplyr::summarise(within_sample_hap = max(.data$within_sample_hap)),
        by = "specimen_name"
      ) |>
      dplyr::rowwise() |>
      dplyr::mutate(within_sample_hap = list(1:.data$within_sample_hap)) |>
      tidyr::unnest("within_sample_hap")

    aa_table_group_only_1_variable_filt_combined <- dplyr::bind_rows(
      aa_table_group_only_1_variable_filt_variable,
      aa_table_group_only_1_variable_filt_invariable
    )

    aa_table_group_only_1_variable_filt_combined_haps <-
      aa_table_group_only_1_variable_filt_combined |>
        dplyr::arrange(
          .data$gene_id, .data$aa_position, .data$within_sample_hap
        ) |>
        dplyr::group_by(
          .data$specimen_name, .data$gene_id, .data$within_sample_hap
        ) |>
        dplyr::summarise(
          positions = paste0(.data$aa_position, collapse = "_"),
          aas = paste0(.data$aa, collapse = "_"),
          wsaf = ifelse(
            all(is.na(.data$wsaf)),
            NA,
            min(.data$wsaf, na.rm = TRUE)
          )
        ) |>
        tidyr::unite(
          "per_gene_variant", "gene_id", "positions", "aas",
          sep = ":"
        ) |>
        dplyr::group_by(.data$specimen_name, .data$within_sample_hap) |>
        dplyr::summarise(
          variant = paste0(.data$per_gene_variant, collapse = ";"),
          wsaf = min(.data$wsaf, na.rm = TRUE)
        )
  } else {
    aa_table_group_only_1_variable_filt_combined_haps <- tibble::tibble()
  }

  aa_table_group_filt <- aa_table_group |>
    dplyr::filter(
      !.data$specimen_name %in% aa_table_group_only_1_variable$specimen_name
    )

  aa_table_group_filt_dominant <- aa_table_group_filt |>
    dplyr::group_by(.data$specimen_name, .data$gene, .data$gene_id, .data$aa_position) |>
    dplyr::mutate(total_reads = sum(.data$reads)) |>
    dplyr::mutate(wsaf = .data$reads / .data$total_reads) |>
    dplyr::filter(.data$wsaf >= wsaf_cut_off) |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::mutate(
      loci_called = dplyr::n_distinct(paste0(.data$gene_id, "-", .data$aa_position))
    ) |>
    dplyr::group_by(.data$specimen_name, .data$gene, .data$gene_id, .data$aa_position) |>
    dplyr::mutate(n_aa = dplyr::n_distinct(.data$aa)) |>
    dplyr::filter(
      all(.data$n_aa == 1),
      .data$loci_called == group_df$loci_in_group[[1]]
    )

  aa_table_group_filt_dominant_collapse <- aa_table_group_filt_dominant |>
    dplyr::arrange(.data$gene_id, .data$aa_position) |>
    dplyr::group_by(.data$specimen_name, .data$gene_id) |>
    dplyr::summarise(
      positions = paste0(.data$aa_position, collapse = "_"),
      aas = paste0(.data$aa, collapse = "_"),
      wsaf = min(.data$wsaf)
    ) |>
    tidyr::unite("per_gene_variant", "gene_id", "positions", "aas", sep = ":") |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::summarise(
      variant = paste0(.data$per_gene_variant, collapse = ";"),
      wsaf = min(.data$wsaf)
    ) |>
    dplyr::group_by(.data$specimen_name) |>
    dplyr::mutate(within_sample_hap = dplyr::row_number())

  dplyr::bind_rows(
    aa_table_group_filt_dominant_collapse,
    aa_table_group_only_1_variable_filt_combined_haps
  )
}

#' Estimate multilocus prevalence and frequency with naive phasing
#'
#' For each loci group, specimens that have calls at every locus are retained.
#' Haplotypes are inferred when exactly one locus is heterozygous, or when every
#' locus has a single allele above `wsaf_cut_off` (including monoclonal
#' samples). Prevalence and frequency are then estimated with `wsaf_prop` or
#' `presence_absence`.
#'
#' @param aa_table Path to a TSV of amino acid calls with columns
#'   `specimen_name`, `gene`, `gene_id`, `aa_position`, `reads`, and `aa`.
#' @param loci_groups_input Path to a TSV of loci groups with columns
#'   `group_id`, `gene_id`, and `aa_position`.
#' @param output_path Optional path for the multilocus prev/freq TSV. If
#'   `NULL`, results are returned without writing.
#' @param recalc_single_locus_output_path Optional path for single-locus
#'   prev/freq recalculated from the inferred multilocus calls.
#' @param method `"wsaf_prop"` (default) or `"presence_absence"`.
#' @param wsaf_cut_off WSAF threshold used to infer a dominant haplotype when
#'   more than one locus is heterozygous. Default: `0.70`.
#'
#' @return A tibble of multilocus prevalence and frequency estimates.
#'
#' @examples
#' aa_path <- system.file(
#'   "extdata", "example2_amino_acid_calls.tsv",
#'   package = "PGEcore"
#' )
#' groups_path <- system.file(
#'   "extdata", "example_loci_groups.tsv",
#'   package = "PGEcore"
#' )
#' multilocus_prevfreq_naive(aa_path, groups_path)
#'
#' @export
multilocus_prevfreq_naive <- function(aa_table,
                                      loci_groups_input,
                                      output_path = NULL,
                                      recalc_single_locus_output_path = NULL,
                                      method = "wsaf_prop",
                                      wsaf_cut_off = 0.70) {
  options(dplyr.summarise.inform = FALSE)
  options(readr.show_col_types = FALSE)

  if (!method %in% c("wsaf_prop", "presence_absence")) {
    stop(method, " is not a valid method", call. = FALSE)
  }
  stopifnot(is.numeric(wsaf_cut_off), length(wsaf_cut_off) == 1L)

  aa_calls <- read_mlp_naive_aa_table(aa_table)
  loci_groups <- read_mlp_naive_loci_groups(loci_groups_input) |>
    dplyr::group_by(.data$group_id) |>
    dplyr::mutate(
      loci_in_group = dplyr::n_distinct(paste0(.data$gene_id, "-", .data$aa_position))
    )

  loci_groups_split <- split(loci_groups, loci_groups$group_id)

  all_calls <- tibble::tibble()
  all_prev_freq <- tibble::tibble()

  for (loci_group in names(loci_groups_split)) {
    group_calls <- build_naive_multilocus_calls_for_group(
      aa_table = aa_calls,
      group_df = loci_groups_split[[loci_group]],
      wsaf_cut_off = wsaf_cut_off
    )

    all_calls <- dplyr::bind_rows(
      all_calls,
      dplyr::mutate(group_calls, group_id = loci_group)
    )

    group_prev_freq <- switch(
      method,
      wsaf_prop = calculate_multilocus_af_prev_wsaf_prop(group_calls),
      presence_absence = calculate_multilocus_af_prev_presence_absence(group_calls)
    ) |>
      dplyr::mutate(group_id = loci_group)

    all_prev_freq <- dplyr::bind_rows(all_prev_freq, group_prev_freq)
  }

  if (!is.null(output_path)) {
    readr::write_tsv(all_prev_freq, output_path)
  }

  if (!is.null(recalc_single_locus_output_path)) {
    slaf_from_ml <- switch(
      method,
      wsaf_prop = generate_single_locus_prev_freq_from_multilocus_groups_wsaf_prop(
        all_calls
      ),
      presence_absence = generate_single_locus_prev_freq_from_multilocus_groups_presence_absence(
        all_calls
      )
    )
    readr::write_tsv(slaf_from_ml, recalc_single_locus_output_path)
  }

  all_prev_freq
}
