#' Add reference sequences extracted from a genome FASTA onto a panel BED table
#'
#' For each BED interval, extracts `length` bases starting at 0-based `start`
#' (Biostrings 1-based `start + 1`) from the contig named in `#chrom`. Minus-
#' strand intervals are reverse-complemented. Requires **Biostrings** (Suggests).
#' Genome FASTA names are truncated at the first whitespace.
#'
#' @param ref_bed Path to a TSV BED table with a header, or a data frame, with
#'   columns `#chrom`, `start`, `end`, `target_name`, `length`, and `strand`.
#' @param genome_fasta Path to a genome FASTA.
#' @param out Optional output TSV path. If `NULL`, results are returned without
#'   writing.
#' @param overwrite If `FALSE` (default), refuse to overwrite `out`.
#'
#' @return A tibble of `ref_bed` with a `ref_seq` column.
#'
#' @export
add_ref_seqs_with_full_genome_ref_fasta <- function(ref_bed,
                                                    genome_fasta,
                                                    out = NULL,
                                                    overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  stop_if_output_exists(out, overwrite)
  check_suggested_pkg(
    "Biostrings",
    "extracting intervals from a genome FASTA via add_ref_seqs_with_full_genome_ref_fasta()"
  )

  if (is.character(ref_bed) && length(ref_bed) == 1L) {
    if (!file.exists(ref_bed)) {
      stop(ref_bed, " does not exist", call. = FALSE)
    }
    bed <- read_ref_bed_table(ref_bed)
    bed_label <- ref_bed
  } else if (is.data.frame(ref_bed)) {
    validate_ref_bed_table(ref_bed, "ref_bed")
    bed <- tibble::as_tibble(ref_bed)
    bed_label <- "ref_bed"
  } else {
    stop("`ref_bed` must be a file path or a data frame.", call. = FALSE)
  }

  stop_on_duplicate_names(bed$target_name, bed_label)

  loaded_genome <- read_genome_dna_string_set(
    genome_fasta,
    "extracting intervals from a genome FASTA via add_ref_seqs_with_full_genome_ref_fasta()"
  )

  bed$ref_seq <- ""
  for (row in seq_len(nrow(bed))) {
    chrom <- bed[["#chrom"]][row]
    if (!chrom %in% names(loaded_genome)) {
      stop(
        chrom, " not in ", genome_fasta, " options: ",
        paste(names(loaded_genome), collapse = ","),
        call. = FALSE
      )
    }
    ref_seq <- Biostrings::subseq(
      loaded_genome[chrom],
      bed$start[row] + 1,
      width = bed$length[row]
    )
    if (identical(bed$strand[row], "-")) {
      ref_seq <- Biostrings::reverseComplement(ref_seq)
    }
    bed$ref_seq[row] <- unname(as.character(ref_seq))
  }

  if (!is.null(out)) {
    readr::write_tsv(bed, out)
  }
  bed
}
