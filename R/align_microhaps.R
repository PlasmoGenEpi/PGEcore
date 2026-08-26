#' Map a gapless sequence coordinate to an aligned coordinate
#'
#' @param aligned_seq Aligned sequence that may contain `"-"` characters.
#' @param pos 1-based position in the gapless sequence.
#' @return The corresponding 1-based index in `aligned_seq`.
#' @keywords internal
get_aln_pos_per_real_pos <- function(aligned_seq, pos) {
  gapless_pos <- which(strsplit(as.character(aligned_seq), NULL)[[1]] != "-")
  gapless_pos[pos]
}

#' Aligned pattern sequence including end gaps from an overlap alignment
#'
#' @param pw_overlapAlign A [pwalign::pairwiseAlignment()] overlap result.
#' @return Character string of the pattern with terminal gaps relative to the
#'   subject.
#' @keywords internal
get_aligned_pattern_from_overlap_align <- function(pw_overlapAlign) {
  overlap_align_aligned_pattern <- pwalign::alignedPattern(pw_overlapAlign)
  overlap_align_aligned_pattern_bases <- as.character(overlap_align_aligned_pattern)

  subject_unaligned <- as.character(pw_overlapAlign@subject@unaligned)
  pattern_unaligned <- as.character(pw_overlapAlign@pattern@unaligned)

  preceding_gap_size <- pw_overlapAlign@subject@range@start - 1
  trailing_gap_size <- nchar(subject_unaligned) -
    (pw_overlapAlign@subject@range@start + pw_overlapAlign@subject@range@width) + 1

  preceding_other_gap_size <- pw_overlapAlign@pattern@range@start - 1
  trailing_other_gap_size <- nchar(pattern_unaligned) -
    (pw_overlapAlign@pattern@range@start + pw_overlapAlign@pattern@range@width) + 1

  if (preceding_other_gap_size > 0) {
    overlap_align_aligned_pattern_bases <- paste0(
      substr(pattern_unaligned, 1, preceding_other_gap_size),
      overlap_align_aligned_pattern_bases
    )
  }
  if (trailing_other_gap_size > 0) {
    overlap_align_aligned_pattern_bases <- paste0(
      overlap_align_aligned_pattern_bases,
      substr(
        pattern_unaligned,
        nchar(pattern_unaligned) - trailing_other_gap_size + 1,
        nchar(pattern_unaligned)
      )
    )
  }
  paste0(
    paste0(rep("-", preceding_gap_size), collapse = ""),
    overlap_align_aligned_pattern_bases,
    paste0(rep("-", trailing_gap_size), collapse = "")
  )
}

#' Aligned subject sequence including end gaps from an overlap alignment
#'
#' @param pw_overlapAlign A [pwalign::pairwiseAlignment()] overlap result.
#' @return Character string of the subject with terminal gaps relative to the
#'   pattern.
#' @keywords internal
get_aligned_subject_from_overlap_align <- function(pw_overlapAlign) {
  overlap_align_aligned_subject <- pwalign::alignedSubject(pw_overlapAlign)
  overlap_align_aligned_subject_bases <- as.character(overlap_align_aligned_subject)

  pattern_unaligned <- as.character(pw_overlapAlign@pattern@unaligned)
  subject_unaligned <- as.character(pw_overlapAlign@subject@unaligned)

  preceding_gap_size <- pw_overlapAlign@pattern@range@start - 1
  trailing_gap_size <- nchar(pattern_unaligned) -
    (pw_overlapAlign@pattern@range@start + pw_overlapAlign@pattern@range@width) + 1

  preceding_other_gap_size <- pw_overlapAlign@subject@range@start - 1
  trailing_other_gap_size <- nchar(subject_unaligned) -
    (pw_overlapAlign@subject@range@start + pw_overlapAlign@subject@range@width) + 1

  if (preceding_other_gap_size > 0) {
    overlap_align_aligned_subject_bases <- paste0(
      substr(subject_unaligned, 1, preceding_other_gap_size),
      overlap_align_aligned_subject_bases
    )
  }
  if (trailing_other_gap_size > 0) {
    overlap_align_aligned_subject_bases <- paste0(
      overlap_align_aligned_subject_bases,
      substr(
        subject_unaligned,
        nchar(subject_unaligned) - trailing_other_gap_size + 1,
        nchar(subject_unaligned)
      )
    )
  }
  paste0(
    paste0(rep("-", preceding_gap_size), collapse = ""),
    overlap_align_aligned_subject_bases,
    paste0(rep("-", trailing_gap_size), collapse = "")
  )
}

#' Require Biostrings and pwalign for microhaplotype overlap alignment
#'
#' @param reason Passed to [check_suggested_pkg()].
#' @return Invisibly `TRUE`.
#' @keywords internal
check_microhap_alignment_pkgs <- function(reason) {
  check_suggested_pkg("Biostrings", reason)
  check_suggested_pkg("pwalign", reason)
  invisible(TRUE)
}

#' Default nucleotide substitution matrix used by pileup and translate
#'
#' @return A substitution matrix from [pwalign::nucleotideSubstitutionMatrix()].
#' @keywords internal
nucleotide_overlap_substitution_matrix <- function() {
  pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -2, baseOnly = TRUE)
}

#' Overlap-align an allele DNAString to a reference sequence
#'
#' @param allele_seq A `DNAString` (already oriented to the plus strand of the
#'   target).
#' @param ref_seq Character reference sequence (already oriented).
#' @param mat Substitution matrix from `nucleotide_overlap_substitution_matrix()`.
#' @return A pairwise alignment object.
#' @keywords internal
overlap_align_allele_to_ref <- function(allele_seq, ref_seq, mat) {
  pwalign::pairwiseAlignment(
    allele_seq,
    Biostrings::DNAString(ref_seq),
    substitutionMatrix = mat,
    gapOpening = 5,
    gapExtension = 1,
    type = "overlap"
  )
}

#' Reverse-complement `ref_seq` for minus-strand targets in a lookup list
#'
#' Operates on a copy so repeated calls do not keep reverse-complementing.
#'
#' @param ref_bed_by_loci_lookup Named list of one-row ref_bed tibbles.
#' @return A copy of the lookup with minus-strand `ref_seq` reverse-complemented.
#' @keywords internal
orient_ref_lookup_to_plus_strand <- function(ref_bed_by_loci_lookup) {
  ref_bed_by_loci_lookup_copy <- ref_bed_by_loci_lookup
  for (target_name in names(ref_bed_by_loci_lookup_copy)) {
    if ("-" == ref_bed_by_loci_lookup_copy[[target_name]]$strand[1]) {
      ref_bed_by_loci_lookup_copy[[target_name]]$ref_seq[1] <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAString(
            ref_bed_by_loci_lookup_copy[[target_name]]$ref_seq[1]
          )
        )
      )
    }
  }
  ref_bed_by_loci_lookup_copy
}

#' Orient an allele sequence to the plus strand of its target
#'
#' @param seq Character allele sequence.
#' @param strand Target strand (`"+"` or `"-"`).
#' @return A `DNAString`.
#' @keywords internal
orient_allele_seq_to_plus_strand <- function(seq, strand) {
  allele_seq <- Biostrings::DNAString(seq)
  if ("-" == strand) {
    allele_seq <- Biostrings::reverseComplement(allele_seq)
  }
  allele_seq
}

#' Index a ref_bed table by `target_name`
#'
#' @param ref_bed Panel location table.
#' @return Named list of one-row tibbles.
#' @keywords internal
ref_bed_lookup_by_target <- function(ref_bed) {
  ref_bed_by_loci <- list()
  for (row in seq_len(nrow(ref_bed))) {
    ref_bed_by_loci[[ref_bed$target_name[row]]] <- ref_bed[row, ]
  }
  ref_bed_by_loci
}

#' Unique haplotypes by target_name and seq
#'
#' @param allele_table Allele table with `target_name` and `seq`.
#' @return Distinct rows sorted by `target_name`.
#' @keywords internal
unique_haps_from_allele_table <- function(allele_table) {
  allele_table |>
    dplyr::select("target_name", "seq") |>
    unique() |>
    dplyr::arrange(.data$target_name)
}

#' Mark which feature rows are fully contained in each panel interval
#'
#' @param ref_bed_tab Panel location table.
#' @param feature_tab Feature BED table (`#chrom`, `start`, `end`).
#' @param out_col Name of the column storing comma-separated feature row numbers.
#' @return `ref_bed_tab` with `out_col` added.
#' @keywords internal
add_intersected_features_to_ref_bed <- function(ref_bed_tab,
                                                feature_tab,
                                                out_col) {
  ref_bed_tab[[out_col]] <- ""
  for (row in seq_len(nrow(ref_bed_tab))) {
    for (feat_row in seq_len(nrow(feature_tab))) {
      if (ref_bed_tab$`#chrom`[row] == feature_tab$`#chrom`[feat_row] &&
          feature_tab$start[feat_row] >= ref_bed_tab$start[row] &&
          feature_tab$end[feat_row] <= ref_bed_tab$end[row]) {
        if ("" != ref_bed_tab[[out_col]][row]) {
          ref_bed_tab[[out_col]][row] <- paste0(ref_bed_tab[[out_col]][row], ",")
        }
        ref_bed_tab[[out_col]][row] <- paste0(
          ref_bed_tab[[out_col]][row],
          feat_row
        )
      }
    }
  }
  ref_bed_tab
}

#' Add `covered_by_target` listing panel targets that fully contain each feature
#'
#' @param feature_tab Feature BED table.
#' @param ref_bed_tab Panel location table.
#' @return `feature_tab` with `covered_by_target` (`"uncovered"` if none).
#' @keywords internal
add_covered_by_target_to_features <- function(feature_tab, ref_bed_tab) {
  feature_tab$covered_by_target <- ""
  for (row in seq_len(nrow(ref_bed_tab))) {
    for (feat_row in seq_len(nrow(feature_tab))) {
      if (ref_bed_tab$`#chrom`[row] == feature_tab$`#chrom`[feat_row] &&
          feature_tab$start[feat_row] >= ref_bed_tab$start[row] &&
          feature_tab$end[feat_row] <= ref_bed_tab$end[row]) {
        if ("" != feature_tab$covered_by_target[feat_row]) {
          feature_tab$covered_by_target[feat_row] <- paste0(
            feature_tab$covered_by_target[feat_row],
            ","
          )
        }
        feature_tab$covered_by_target[feat_row] <- paste0(
          feature_tab$covered_by_target[feat_row],
          ref_bed_tab$target_name[row]
        )
      }
    }
  }
  feature_tab |>
    dplyr::mutate(
      covered_by_target = ifelse(
        "" == .data$covered_by_target,
        "uncovered",
        .data$covered_by_target
      )
    )
}

#' Features intersecting a target, with coordinates relative to the target start
#'
#' @param lookup_row One-row ref_bed tibble including the intersect column.
#' @param feature_tab Full feature table.
#' @param intersect_col Column of comma-separated feature row numbers.
#' @return Subset of `feature_tab` with `rel_start` and `rel_end`.
#' @keywords internal
features_for_target_with_rel_coords <- function(lookup_row,
                                                feature_tab,
                                                intersect_col) {
  feature_tab[
    as.numeric(unlist(strsplit(lookup_row[[intersect_col]][1], ","))),
  ] |>
    dplyr::mutate(
      rel_start = .data$start - lookup_row$start[1],
      rel_end = .data$end - lookup_row$start[1]
    )
}

#' Warnings when requested specimen/target names are missing from an allele table
#'
#' @param allele_data Allele table.
#' @param select_target_names Requested targets (empty = none).
#' @param select_specimen_names Requested specimens (empty = none).
#' @param allele_table_fnp Path label used in messages.
#' @return Character vector of warning strings (possibly empty).
#' @keywords internal
warnings_for_subselecting_allele_table <- function(allele_data,
                                                   select_target_names,
                                                   select_specimen_names,
                                                   allele_table_fnp) {
  warns <- character(0)
  if (length(select_specimen_names) > 0) {
    missing_sel_specs <- setdiff(
      select_specimen_names,
      unique(allele_data$specimen_name)
    )
    if (length(missing_sel_specs) > 0) {
      warns <- c(
        warns,
        paste0(
          "supplied --select_specimen_names but the following specimen_names ",
          "are missing from ", allele_table_fnp, "\n",
          paste(missing_sel_specs, collapse = ",")
        )
      )
    }
  }
  if (length(select_target_names) > 0) {
    missing_sel_tars <- setdiff(
      select_target_names,
      unique(allele_data$target_name)
    )
    if (length(missing_sel_tars) > 0) {
      warns <- c(
        warns,
        paste0(
          "supplied --select_target_names but the following target_names ",
          "are missing from ", allele_table_fnp, "\n",
          paste(missing_sel_tars, collapse = ",")
        )
      )
    }
  }
  warns
}

#' Coerce ref_bed / allele_table / feature table inputs to tibbles
#'
#' @param x A path or data frame.
#' @param reader Function used when `x` is a path.
#' @param what Label for error messages.
#' @return A tibble.
#' @keywords internal
as_input_tibble <- function(x, reader, what) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop(x, " does not exist", call. = FALSE)
    }
    return(reader(x))
  }
  if (is.data.frame(x)) {
    return(tibble::as_tibble(x))
  }
  stop("`", what, "` must be a file path or a data frame.", call. = FALSE)
}
