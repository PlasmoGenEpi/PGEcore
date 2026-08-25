#' Extract SNP bases from unique haplotypes via overlap alignment
#'
#' @param allele_table_unique_haps_tab Unique `target_name`/`seq` rows.
#' @param microhaps_intersected_with_snps_of_interest Targets covering SNPs.
#' @param ref_bed_by_loci_lookup Named list of one-row ref_bed tibbles.
#' @param snps_of_interest_tab SNP-of-interest table.
#' @return Tibble of SNP calls per haplotype.
#' @keywords internal
extract_snps_of_interest <- function(allele_table_unique_haps_tab,
                                     microhaps_intersected_with_snps_of_interest,
                                     ref_bed_by_loci_lookup,
                                     snps_of_interest_tab) {
  mat <- nucleotide_overlap_substitution_matrix()
  all_snps <- tibble::tibble(
    target_name = character(),
    seq = character(),
    chrom = character(),
    pos = numeric(),
    snp_name = character(),
    strand = character(),
    ref_base = character(),
    seq_base = character()
  )
  ref_lookup <- orient_ref_lookup_to_plus_strand(ref_bed_by_loci_lookup)

  for (row in seq_len(nrow(allele_table_unique_haps_tab))) {
    target <- allele_table_unique_haps_tab$target_name[row]
    if (!(target %in% microhaps_intersected_with_snps_of_interest)) {
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
    snps_for_target <- features_for_target_with_rel_coords(
      ref_lookup[[target]],
      snps_of_interest_tab,
      "intersected_snps_of_interest"
    )
    snps_for_hap <- tibble::tibble()
    for (snp_row in seq_len(nrow(snps_for_target))) {
      aln_pos <- get_aln_pos_per_real_pos(
        get_aligned_subject_from_overlap_align(overlap_align),
        snps_for_target$rel_start[snp_row] + 1
      )
      aligned_pattern <- get_aligned_pattern_from_overlap_align(overlap_align)
      seq_base <- Biostrings::DNAString(
        substr(
          aligned_pattern,
          aln_pos,
          aln_pos + snps_for_target$length[snp_row] - 1
        )
      )
      ref_base <- Biostrings::DNAString(
        substr(
          ref_lookup[[target]]$ref_seq[1],
          snps_for_target$rel_start[snp_row] + 1,
          snps_for_target$rel_start[snp_row] + 1 +
            snps_for_target$length[snp_row] - 1
        )
      )
      if ("-" == snps_for_target$strand[snp_row]) {
        seq_base <- Biostrings::reverseComplement(seq_base)
        ref_base <- Biostrings::reverseComplement(ref_base)
      }
      snps_for_hap <- dplyr::bind_rows(
        snps_for_hap,
        tibble::tibble(
          target_name = target,
          seq = allele_table_unique_haps_tab$seq[row],
          chrom = snps_for_target$`#chrom`[snp_row],
          pos = snps_for_target$start[snp_row],
          snp_name = snps_for_target$name[snp_row],
          strand = snps_for_target$strand[snp_row],
          ref_base = as.character(ref_base),
          seq_base = as.character(seq_base)
        )
      )
    }
    all_snps <- dplyr::bind_rows(all_snps, snps_for_hap)
  }
  all_snps
}

#' Collapse overlapping-target SNP calls by summing or picking the best target
#'
#' @param allele_table_to_collapse SNP calls joined to allele counts.
#' @param collapse_calls_by_summing If `TRUE`, sum reads across targets.
#' @return Collapsed tibble.
#' @keywords internal
collapse_snp_allele_table <- function(allele_table_to_collapse,
                                      collapse_calls_by_summing = FALSE) {
  if (isTRUE(collapse_calls_by_summing)) {
    allele_table_to_collapse |>
      dplyr::group_by(
        .data$specimen_name,
        .data$chrom,
        .data$pos,
        .data$snp_name,
        .data$strand,
        .data$ref_base,
        .data$seq_base
      ) |>
      dplyr::summarise(
        reads = sum(.data$reads),
        target_name = paste0(sort(.data$target_name), collapse = ",")
      )
  } else {
    winner <- allele_table_to_collapse |>
      dplyr::group_by(
        .data$specimen_name,
        .data$chrom,
        .data$pos,
        .data$snp_name,
        .data$ref_base,
        .data$target_name
      ) |>
      dplyr::summarise(reads = sum(.data$reads)) |>
      dplyr::arrange(dplyr::desc(.data$reads)) |>
      dplyr::mutate(
        reads_rank = dplyr::row_number(),
        covered_by_target_names = paste0(sort(.data$target_name), collapse = ",")
      ) |>
      dplyr::filter(.data$reads_rank == 1) |>
      dplyr::ungroup() |>
      dplyr::select(-"reads_rank") |>
      dplyr::rename(best_target_name = "target_name")

    collapsed <- allele_table_to_collapse |>
      dplyr::left_join(
        winner |> dplyr::ungroup() |> dplyr::select(-"reads"),
        by = c("specimen_name", "chrom", "pos", "snp_name", "ref_base")
      ) |>
      dplyr::filter(.data$target_name == .data$best_target_name) |>
      dplyr::select(-"seq")

    collapsed |>
      dplyr::group_by(
        .data$specimen_name,
        .data$target_name,
        .data$chrom,
        .data$pos,
        .data$snp_name,
        .data$strand,
        .data$ref_base,
        .data$seq_base,
        .data$best_target_name,
        .data$covered_by_target_names
      ) |>
      dplyr::summarise(reads = sum(.data$reads))
  }
}

#' Pile up specific SNPs covered by microhaplotype sequences
#'
#' Aligns each unique haplotype to its panel reference with an overlap
#' pairwise alignment, then extracts bases at SNP-of-interest coordinates.
#'
#' ## Inputs
#'
#' - **`allele_table`**: Allele table (`specimen_name`, `target_name`, `reads`,
#'   `seq`), as a file path or data frame. See
#'   `vignette("input-formats", package = "PGEcore")`.
#' - **`ref_bed`**: Panel BED with `ref_seq` (`#chrom`, `start`, `end`,
#'   `target_name`, `length`, `strand`, `ref_seq`).
#' - **`snps_of_interest`**: SNP BED (`#chrom`, `start`, `end`, `name`,
#'   `length`, `strand`); each SNP must span one base (`end - start == 1`).
#'
#' ## Outputs
#'
#' - **`output_dir`**: Directory receiving `snp_calls.tsv.gz`,
#'   `collapsed_snp_calls.tsv.gz`, `snps_covered_by_target_samples_info.tsv`,
#'   and optionally `allele_table_out_uncallable.tsv`.
#'
#' ## Running
#'
#' ```r
#' pileup_specific_snps(
#'   allele_table = "allele_table.tsv",
#'   ref_bed = "ref_bed_with_seq.tsv",
#'   snps_of_interest = "snps.bed",
#'   output_dir = "pileup_out"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/pileup_specific_snps \
#'   --allele_table allele_table.tsv \
#'   --ref_bed ref_bed_with_seq.tsv \
#'   --snps_of_interest snps.bed \
#'   --output_dir pileup_out
#' ```
#'
#' Requires **Biostrings** and **pwalign** (Suggests).
#'
#' @param allele_table Path or data frame of allele table. See *Inputs*.
#' @param ref_bed Path or data frame of panel BED with `ref_seq`. See *Inputs*.
#' @param snps_of_interest Path or data frame of SNP BED. See *Inputs*.
#' @param output_dir Directory to write results. Created if missing.
#' @param select_target_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of targets to keep.
#' @param select_specimen_names Optional comma-separated names, path to a
#'   one-column TSV, or character vector of specimens to keep.
#' @param overwrite_dir If `FALSE` (default), refuse to replace an existing
#'   `output_dir`.
#' @param collapse_calls_by_summing If `TRUE`, sum reads across overlapping
#'   targets; otherwise keep the target with the highest read count.
#'
#' @return A named list with `snp_calls`, `collapsed_snp_calls`,
#'   `snps_covered_by_target_samples_info`, and `uncallable`.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
pileup_specific_snps <- function(allele_table,
                                 ref_bed,
                                 snps_of_interest,
                                 output_dir,
                                 select_target_names = NULL,
                                 select_specimen_names = NULL,
                                 overwrite_dir = FALSE,
                                 collapse_calls_by_summing = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  options(readr.show_col_types = FALSE)
  check_microhap_alignment_pkgs("pileup_specific_snps()")
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
  if (is.character(snps_of_interest) && length(snps_of_interest) == 1L) {
    snp_label <- snps_of_interest
    snps_of_interest <- as_input_tibble(
      snps_of_interest,
      function(path) {
        tab <- readr::read_tsv(path, col_names = TRUE, show_col_types = FALSE)
        validate_required_columns(
          tab,
          c("#chrom", "start", "end", "name", "length", "strand"),
          path
        )
        tab
      },
      "snps_of_interest"
    )
  } else {
    snp_label <- "snps_of_interest"
    snps_of_interest <- as_input_tibble(
      snps_of_interest,
      identity,
      "snps_of_interest"
    )
    validate_required_columns(
      snps_of_interest,
      c("#chrom", "start", "end", "name", "length", "strand"),
      snp_label
    )
  }
  validate_ref_bed_with_seq_table(ref_bed, ref_label)

  warnings <- character(0)
  for (row in seq_len(nrow(snps_of_interest))) {
    if ((snps_of_interest$end[row] - snps_of_interest$start[row]) != 1) {
      warnings <- c(
        warnings,
        paste0(
          "snps of interest must be of length 1, locus: ",
          snps_of_interest$name[row],
          " is length: ",
          snps_of_interest$length[row]
        )
      )
    }
  }
  warnings <- c(
    warnings,
    warnings_for_subselecting_allele_table(
      allele_table,
      select_targets,
      select_specs,
      allele_label
    )
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
        "the following snps were missing from the reference location file ",
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
    snps_of_interest,
    "intersected_snps_of_interest"
  )
  snps_of_interest <- add_covered_by_target_to_features(snps_of_interest, ref_bed)
  ref_bed_by_loci <- ref_bed_lookup_by_target(ref_bed)
  microhaps_with_snps <- ref_bed |>
    dplyr::filter("" != .data$intersected_snps_of_interest) |>
    dplyr::pull("target_name")

  hap_snps <- extract_snps_of_interest(
    allele_table_unique_haps,
    microhaps_with_snps,
    ref_bed_by_loci,
    snps_of_interest
  )

  allele_table_out <- allele_table |>
    dplyr::filter(.data$target_name %in% microhaps_with_snps) |>
    dplyr::left_join(
      hap_snps,
      relationship = "many-to-many",
      by = c("target_name", "seq")
    )

  covered_by_samples <- allele_table_out |>
    dplyr::group_by(.data$chrom, .data$pos, .data$snp_name, .data$ref_base) |>
    dplyr::summarise(
      n_samples = dplyr::n_distinct(.data$specimen_name),
      .groups = "drop"
    ) |>
    dplyr::mutate(total_samples = dplyr::n_distinct(allele_table$specimen_name))

  snps_of_interest_out <- snps_of_interest |>
    dplyr::left_join(
      covered_by_samples |> dplyr::rename(name = "snp_name"),
      by = c("name")
    )

  uncallable <- allele_table_out |>
    dplyr::filter(.data$seq_base == "-")
  allele_table_out_filt <- allele_table_out |>
    dplyr::filter(.data$seq_base != "-")

  collapsed <- collapse_snp_allele_table(
    allele_table_out_filt,
    collapse_calls_by_summing
  )

  allele_table_out <- allele_table_out |>
    dplyr::group_by(.data$chrom, .data$pos, .data$snp_name, .data$ref_base) |>
    dplyr::mutate(is_biallelic = dplyr::n_distinct(.data$seq_base) <= 2)

  collapsed <- collapsed |>
    dplyr::group_by(.data$chrom, .data$pos, .data$snp_name, .data$ref_base) |>
    dplyr::mutate(is_biallelic = dplyr::n_distinct(.data$seq_base) <= 2)

  readr::write_tsv(
    allele_table_out,
    file.path(output_dir, "snp_calls.tsv.gz")
  )
  readr::write_tsv(
    collapsed,
    file.path(output_dir, "collapsed_snp_calls.tsv.gz")
  )
  readr::write_tsv(
    snps_of_interest_out,
    file.path(output_dir, "snps_covered_by_target_samples_info.tsv")
  )
  if (nrow(uncallable) > 0) {
    readr::write_tsv(
      uncallable,
      file.path(output_dir, "allele_table_out_uncallable.tsv")
    )
  }

  list(
    snp_calls = allele_table_out,
    collapsed_snp_calls = collapsed,
    snps_covered_by_target_samples_info = snps_of_interest_out,
    uncallable = uncallable
  )
}
