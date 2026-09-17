#' Add reference sequences extracted from a genome FASTA onto a panel BED table
#'
#' For each BED interval, extracts `length` bases starting at 0-based `start`
#' (Biostrings 1-based `start + 1`) from the contig named in `#chrom`. Minus-
#' strand intervals are reverse-complemented. Genome FASTA names are truncated
#' at the first whitespace.
#'
#' ## Inputs
#'
#' - **`ref_bed`**: Panel BED TSV with header (`#chrom`, `start`, `end`,
#'   `target_name`, `length`, `strand`), as a file path or data frame.
#' - **`genome_fasta`**: Genome FASTA to extract intervals from.
#'
#' ## Outputs
#'
#' - **`output`** (optional): `ref_bed` TSV with a `ref_seq` column. If `NULL`,
#'   results are returned without writing.
#'
#' ## Running
#'
#' ```r
#' add_ref_seqs_with_full_genome_ref_fasta(
#'   ref_bed = "ref_bed.tsv",
#'   genome_fasta = "genome.fasta",
#'   output = "ref_bed_with_seq.tsv"
#' )
#' ```
#'
#' ```bash
#' Rscript exec/add_ref_seqs_with_full_genome_ref_fasta \
#'   --ref_bed ref_bed.tsv \
#'   --genome_fasta genome.fasta \
#'   --output ref_bed_with_seq.tsv
#' ```
#'
#' Requires **Biostrings** (Suggests).
#'
#' @param ref_bed Path to a panel BED TSV, or a data frame with the same
#'   columns. See *Inputs*.
#' @param genome_fasta Path to a genome FASTA. See *Inputs*.
#' @param output Optional output TSV path.
#' @param overwrite If `FALSE` (default), refuse to overwrite `output`.
#'
#' @return A tibble of `ref_bed` with a `ref_seq` column.
#'
#' @seealso `vignette("input-formats", package = "PGEcore")`
#'
#' @export
add_ref_seqs_with_full_genome_ref_fasta <- function(ref_bed,
                                                    genome_fasta,
                                                    output = NULL,
                                                    overwrite = FALSE) {
  options(dplyr.summarise.inform = FALSE)
  stop_if_output_exists(output, overwrite)
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

  if (!is.null(output)) {
    readr::write_tsv(bed, output)
  }
  bed
}
