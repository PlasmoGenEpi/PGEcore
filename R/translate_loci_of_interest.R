#' Translate amino-acid loci from unique haplotypes via overlap alignment
#'
#' @param allele_table_unique_haps_tab Unique `target_name`/`seq` rows.
#' @param microhaps_intersected_with_loci_of_interest Targets covering loci.
#' @param ref_bed_by_loci_lookup Named list of one-row ref_bed tibbles.
#' @param loci_of_interest_tab Loci-of-interest table.
#' @return Tibble of codon/AA calls per haplotype.
#' @keywords internal
translate_microhap_seqs <- function(allele_table_unique_haps_tab,
                                    microhaps_intersected_with_loci_of_interest,
                                    ref_bed_by_loci_lookup,
                                    loci_of_interest_tab) {
  all_loci <- tibble::tibble(
    target_name = character(),
    seq = character(),
    gene = character(),
    gene_id = character(),
    aa_position = numeric(),
    ref_codon = character(),
    ref_aa = character(),
    codon = character(),
    aa = character()
  )
  mat <- nucleotide_overlap_substitution_matrix()
  ref_lookup <- orient_ref_lookup_to_plus_strand(ref_bed_by_loci_lookup)

  for (row in seq_len(nrow(allele_table_unique_haps_tab))) {
    target <- allele_table_unique_haps_tab$target_name[row]
    if (!(target %in% microhaps_intersected_with_loci_of_interest)) {
      next
    }
    allele_seq <- orient_allele_seq_to_plus_strand(
      allele_table_unique_haps_tab$seq[row],
      ref_lookup[[target]]$strand
    )
    overlap_align <- overlap_align_allele_to_ref(
      allele_seq,
      ref_lookup[[target]]$ref_seq[1],
      mat
    )
    loci_for_target <- features_for_target_with_rel_coords(
      ref_lookup[[target]],
      loci_of_interest_tab,
      "intersected_loci_of_interest"
    )
    loci_for_hap <- tibble::tibble()
    for (loci_row in seq_len(nrow(loci_for_target))) {
      aln_pos <- get_aln_pos_per_real_pos(
        get_aligned_subject_from_overlap_align(overlap_align),
        loci_for_target$rel_start[loci_row] + 1
      )
      aligned_pattern <- get_aligned_pattern_from_overlap_align(overlap_align)
      seq_codon <- Biostrings::DNAString(substr(aligned_pattern, aln_pos, aln_pos + 2))
      ref_codon <- Biostrings::DNAString(
        substr(
          ref_lookup[[target]]$ref_seq[1],
          loci_for_target$rel_start[loci_row] + 1,
          loci_for_target$rel_start[loci_row] + 2 + 1
        )
      )
      if ("-" == loci_for_target$strand[loci_row]) {
        seq_codon <- Biostrings::reverseComplement(seq_codon)
        ref_codon <- Biostrings::reverseComplement(ref_codon)
      }
      if (grepl("-", as.character(seq_codon), fixed = TRUE)) {
        seq_aa <- "X"
      } else {
        seq_aa <- Biostrings::translate(seq_codon, no.init.codon = TRUE)
      }
      ref_aa <- Biostrings::translate(ref_codon, no.init.codon = TRUE)
      loci_for_hap <- dplyr::bind_rows(
        loci_for_hap,
        tibble::tibble(
          target_name = target,
          seq = allele_table_unique_haps_tab$seq[row],
          gene = loci_for_target$gene[loci_row],
          gene_id = loci_for_target$gene_id[loci_row],
          aa_position = loci_for_target$aa_position[loci_row],
          ref_codon = as.character(ref_codon),
          ref_aa = as.character(ref_aa),
          codon = as.character(seq_codon),
          aa = as.character(seq_aa)
        )
      )
    }
    all_loci <- dplyr::bind_rows(all_loci, loci_for_hap)
  }
  all_loci |>
    dplyr::mutate(aa_locus = paste0(.data$gene_id, ":", .data$aa_position))
}

#' Collapse overlapping-target amino-acid calls
#'
#' @param allele_table_to_filter Translated calls joined to allele counts.
#' @param collapse_calls_by_summing If `TRUE`, sum reads across targets.
#' @return Collapsed tibble.
#' @keywords internal
collapse_aa_allele_table <- function(allele_table_to_filter,
                                     collapse_calls_by_summing = FALSE) {
  if (isTRUE(collapse_calls_by_summing)) {
    allele_table_to_filter |>
      dplyr::group_by(
        .data$specimen_name,
        .data$gene,
        .data$gene_id,
        .data$aa_position,
        .data$aa_locus,
        .data$ref_aa,
        .data$aa
      ) |>
      dplyr::summarise(
        reads = sum(.data$reads),
        target_name = paste0(unique(sort(.data$target_name)), collapse = ",")
      )
  } else {
    winner <- allele_table_to_filter |>
      dplyr::group_by(
        .data$specimen_name,
        .data$gene,
        .data$gene_id,
        .data$aa_position,
        .data$aa_locus,
        .data$ref_aa,
        .data$target_name
      ) |>
      dplyr::summarise(reads = sum(.data$reads)) |>
      dplyr::arrange(dplyr::desc(.data$reads)) |>
      dplyr::mutate(
        reads_rank = dplyr::row_number(),
        covered_by_target_names = paste0(
          unique(sort(.data$target_name)),
          collapse = ","
        )
      ) |>
      dplyr::filter(.data$reads_rank == 1) |>
      dplyr::ungroup() |>
      dplyr::select(-"reads_rank") |>
      dplyr::rename(best_target_name = "target_name")

    collapsed <- allele_table_to_filter |>
      dplyr::left_join(
        winner |> dplyr::ungroup() |> dplyr::select(-"reads"),
        by = c(
          "specimen_name", "gene", "gene_id", "aa_position", "aa_locus", "ref_aa"
        )
      ) |>
      dplyr::filter(.data$target_name == .data$best_target_name) |>
      dplyr::select(-"seq")

    collapsed |>
      dplyr::group_by(
        .data$specimen_name,
        .data$target_name,
        .data$gene,
        .data$gene_id,
        .data$aa_position,
        .data$aa_locus,
        .data$ref_aa,
        .data$aa,
        .data$best_target_name,
        .data$covered_by_target_names
      ) |>
      dplyr::summarise(reads = sum(.data$reads))
  }
}

#' Validate column types for translate_loci_of_interest inputs
#'
#' @param ref_bed Panel BED with `ref_seq`.
#' @param loci_of_interest Loci-of-interest table.
#' @param allele_table Allele table.
#' @return Character vector of warning strings (possibly empty).
#' @keywords internal
validate_translate_column_types <- function(ref_bed,
                                            loci_of_interest,
                                            allele_table) {
  warns <- character(0)
  ref_bed_rules <- validate::validator(
    is.character(`#chrom`),
    is.numeric(start),
    is.numeric(end),
    is.character(target_name),
    is.numeric(length),
    is.character(strand),
    is.character(ref_seq),
    !is.na(`#chrom`),
    !is.na(start),
    !is.na(end),
    !is.na(target_name),
    !is.na(length),
    !is.na(strand),
    !is.na(ref_seq)
  )
  msg <- warn_on_validate_fails(ref_bed, ref_bed_rules, "ref_bed")
  if (!is.null(msg)) {
    warns <- c(warns, msg)
  }

  loci_rules <- validate::validator(
    is.character(`#chrom`),
    is.numeric(start),
    is.numeric(end),
    is.character(name),
    is.numeric(length),
    is.character(strand),
    is.character(gene),
    is.numeric(aa_position),
    is.character(gene_id),
    !is.na(`#chrom`),
    !is.na(start),
    !is.na(end),
    !is.na(name),
    !is.na(length),
    !is.na(strand),
    !is.na(gene),
    !is.na(aa_position),
    !is.na(gene_id)
  )
  msg <- warn_on_validate_fails(loci_of_interest, loci_rules, "loci_of_interest")
  if (!is.null(msg)) {
    warns <- c(warns, msg)
  }

  allele_rules <- validate::validator(
    is.character(specimen_name),
    is.numeric(reads),
    is.character(target_name),
    is.character(seq),
    !is.na(specimen_name),
    !is.na(reads),
    !is.na(target_name),
    !is.na(seq)
  )
  msg <- warn_on_validate_fails(allele_table, allele_rules, "allele_table")
  if (!is.null(msg)) {
    warns <- c(warns, msg)
  }
  warns
}

#' Translate loci of interest from microhaplotype sequences
#'
#' Aligns each unique haplotype to its panel reference with an overlap
#' pairwise alignment, extracts the codon at each locus of interest, and
#' translates it.
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele table (`specimen_name`, `target_name`, `reads`,
#'   `seq`), as a file path or data frame. See
#'   `vignette("input-formats", package = "PGEcore")`.
#' - **`ref_bed`**: Panel BED with `ref_seq` (`#chrom`, `start`, `end`,
#'   `target_name`, `length`, `strand`, `ref_seq`).
#' - **`loci_of_interest`**: Codon BED (`#chrom`, `start`, `end`, `name`,
#'   `length`, `strand`, `gene`, `gene_id`, `aa_position`); each locus must have
#'   `length == 3`.
#'
#' ## Outputs
#'
#' - **`output_dir`**: Directory receiving
#'   `loci_of_interest_for_target_for_microhap.tsv.gz`,
#'   `amino_acid_calls.tsv.gz`, `collapsed_amino_acid_calls.tsv.gz`,
#'   `loci_covered_by_target_samples_info.tsv`, and optionally
#'   `allele_table_out_untranslatable.tsv`.
#'
#' ## Running
#'
#' ```r
#' translate_loci_of_interest(
#'   allele_table = "allele_table.tsv",
#'   ref_bed = "ref_bed_with_seq.tsv",
#'   loci_of_interest = "loci.bed",
#'   output_dir = "translate_out"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/translate_loci_of_interest \
#'   --allele_table allele_table.tsv \
#'   --ref_bed ref_bed_with_seq.tsv \
#'   --loci_of_interest loci.bed \
#'   --output_dir translate_out
#' ```
#'
#' Requires **Biostrings** and **pwalign** (Suggests).
#'
#' @param allele_table Path or data frame of allele table. See *Inputs*.
#' @param ref_bed Path or data frame of panel BED with `ref_seq`. See *Inputs*.
#' @param loci_of_interest Path or data frame of codon BED. See *Inputs*.
#' @param output_dir Directory to write results. Created if missing.
#' @param select_target_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of targets to keep.
#' @param select_specimen_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of specimens to keep.
#' @param overwrite_dir If `FALSE` (default), refuse to replace an existing
#'   `output_dir`.
#' @param output_stop_codons If `FALSE` (default), treat `*` as untranslatable
#'   (along with `X`).
#' @param collapse_calls_by_summing If `TRUE`, sum reads across overlapping
#'   targets; otherwise keep the target with the highest read count.
#'
#' @return A named list with `loci_of_interest_for_target_for_microhap`,
#'   `amino_acid_calls`, `collapsed_amino_acid_calls`,
#'   `loci_covered_by_target_samples_info`, and `untranslatable`.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
translate_loci_of_interest <- function(allele_table,
                                       ref_bed,
                                       loci_of_interest,
                                       output_dir,
                                       select_target_names = NULL,
                                       select_specimen_names = NULL,
                                       overwrite_dir = FALSE,
                                       output_stop_codons = FALSE,
                                       collapse_calls_by_summing = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  options(readr.show_col_types = FALSE)
  check_microhap_alignment_pkgs("translate_loci_of_interest()")
  ensure_output_directory(output_dir, overwrite_dir)

  select_targets <- parse_name_list_arg(select_target_names)
  select_specs <- parse_name_list_arg(select_specimen_names)

  allele_label <- if (is.character(allele_table) && length(allele_table) == 1L) {
    allele_table
  } else {
    "allele_table"
  }
  ref_label <- if (is.character(ref_bed) && length(ref_bed) == 1L) {
    ref_bed
  } else {
    "ref_bed"
  }

  allele_table <- as_input_tibble(
    allele_table,
    read_mhap_allele_table,
    "allele_table"
  )
  validate_required_columns(
    allele_table,
    c("specimen_name", "target_name", "reads", "seq"),
    allele_label
  )
  ref_bed <- as_input_tibble(ref_bed, read_ref_bed_with_seq_table, "ref_bed")
  loci_reader <- function(path) {
    tab <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
    validate_required_columns(
      tab,
      c(
        "#chrom", "start", "end", "name", "length", "strand",
        "gene", "gene_id", "aa_position"
      ),
      path
    )
    tab
  }
  if (is.character(loci_of_interest) && length(loci_of_interest) == 1L) {
    loci_of_interest <- as_input_tibble(
      loci_of_interest,
      loci_reader,
      "loci_of_interest"
    )
  } else {
    loci_of_interest <- as_input_tibble(
      loci_of_interest,
      identity,
      "loci_of_interest"
    )
    validate_required_columns(
      loci_of_interest,
      c(
        "#chrom", "start", "end", "name", "length", "strand",
        "gene", "gene_id", "aa_position"
      ),
      "loci_of_interest"
    )
  }
  validate_ref_bed_with_seq_table(ref_bed, ref_label)

  warnings <- warnings_for_subselecting_allele_table(
    allele_table,
    select_targets,
    select_specs,
    allele_label
  )
  warnings <- c(
    warnings,
    validate_translate_column_types(ref_bed, loci_of_interest, allele_table)
  )

  if (length(select_targets) > 0) {
    missing_sel_tars <- setdiff(select_targets, ref_bed$target_name)
    if (length(missing_sel_tars) > 0) {
      warnings <- c(
        warnings,
        paste0(
          "supplied --select_target_names but the following targets are ",
          "missing from ", ref_label, "\n",
          paste(missing_sel_tars, collapse = ",")
        )
      )
    }
    ref_bed <- dplyr::filter(ref_bed, .data$target_name %in% select_targets)
  }

  for (row in seq_len(nrow(loci_of_interest))) {
    if (loci_of_interest$length[row] != 3) {
      warnings <- c(
        warnings,
        paste0(
          "loci of interest must be of length 3, locus: ",
          loci_of_interest$name[row],
          " is length: ",
          loci_of_interest$length[row]
        )
      )
    }
  }

  allele_table <- filter_snp_table_for_optional_subselecting(
    allele_table,
    select_targets,
    select_specs
  )
  ref_allele_decomp <- set_decompose(
    ref_bed$target_name,
    unique(allele_table$target_name)
  )
  if (length(ref_allele_decomp$only_in_vector_b) > 0) {
    warnings <- c(
      warnings,
      paste0(
        "the following loci were missing from the reference location file ",
        ref_label, " but are in ", allele_label, "\n",
        paste(ref_allele_decomp$only_in_vector_b, collapse = ",")
      )
    )
  }
  if (length(warnings) > 0) {
    stop(paste0("\n", paste(warnings, collapse = "\n")), call. = FALSE)
  }

  allele_table_unique_haps <- unique_haps_from_allele_table(allele_table)
  ref_bed <- add_intersected_features_to_ref_bed(
    ref_bed,
    loci_of_interest,
    "intersected_loci_of_interest"
  )
  loci_of_interest <- add_covered_by_target_to_features(
    loci_of_interest,
    ref_bed
  )
  microhaps_with_loci <- ref_bed |>
    dplyr::filter("" != .data$intersected_loci_of_interest) |>
    dplyr::pull("target_name")
  ref_bed_by_loci <- ref_bed_lookup_by_target(ref_bed)

  hap_loci <- translate_microhap_seqs(
    allele_table_unique_haps,
    microhaps_with_loci,
    ref_bed_by_loci,
    loci_of_interest
  )

  allele_table_out <- allele_table |>
    dplyr::filter(.data$target_name %in% microhaps_with_loci) |>
    dplyr::left_join(
      hap_loci,
      relationship = "many-to-many",
      by = c("target_name", "seq")
    )

  covered_by_samples <- allele_table_out |>
    dplyr::group_by(
      .data$gene,
      .data$gene_id,
      .data$aa_position,
      .data$ref_aa
    ) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(.data$specimen_name),
      .groups = "drop"
    ) |>
    dplyr::mutate(total_samples = dplyr::n_distinct(allele_table$specimen_name))

  loci_of_interest_out <- loci_of_interest |>
    dplyr::left_join(
      covered_by_samples,
      by = c("gene", "aa_position", "gene_id")
    )

  untranslateable_calls <- c("X")
  if (!isTRUE(output_stop_codons)) {
    untranslateable_calls <- c(untranslateable_calls, "*")
  }
  untranslatable <- allele_table_out |>
    dplyr::filter(.data$aa %in% untranslateable_calls)
  allele_table_out_filt <- allele_table_out |>
    dplyr::filter(!(.data$aa %in% untranslateable_calls))

  collapsed <- collapse_aa_allele_table(
    allele_table_out_filt,
    collapse_calls_by_summing
  )

  readr::write_tsv(
    hap_loci,
    file.path(output_dir, "loci_of_interest_for_target_for_microhap.tsv.gz")
  )
  readr::write_tsv(
    allele_table_out,
    file.path(output_dir, "amino_acid_calls.tsv.gz")
  )
  readr::write_tsv(
    collapsed,
    file.path(output_dir, "collapsed_amino_acid_calls.tsv.gz")
  )
  readr::write_tsv(
    loci_of_interest_out,
    file.path(output_dir, "loci_covered_by_target_samples_info.tsv")
  )
  if (nrow(untranslatable) > 0) {
    readr::write_tsv(
      untranslatable,
      file.path(output_dir, "allele_table_out_untranslatable.tsv")
    )
  }

  list(
    loci_of_interest_for_target_for_microhap = hap_loci,
    amino_acid_calls = allele_table_out,
    collapsed_amino_acid_calls = collapsed,
    loci_covered_by_target_samples_info = loci_of_interest_out,
    untranslatable = untranslatable
  )
}
